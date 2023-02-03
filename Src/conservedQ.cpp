#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(false);

    std::string plotFileName; pp.get("infile",plotFileName);
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    //init_mech();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Vector<std::string> varNames = {"density","z_velocity","rhoh","temp","Y(H2)"};

    const int nCompIn  = varNames.size();
    int nCompOut = 13;
    
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
    }
    //in
    const int idrholocal = 0; //rho - keep
    const int iduzlocal = 1; //zvel - keep
    const int idrhohlocal = 2; //rhoh - keep
    const int idtemplocal = 3; //temp - keep
    const int idh2local = 4; //h2
    //out
    const int idrhouzlocal = 4; // rho*uz
    const int idrhouzuzlocal = 5; // rho*uz*uz
    const int idrhohuzlocal = 6; //rho*h*uz
    const int idrhohhlocal = 7; //rho*h*h
    const int idtempuzlocal = 8; // T*uz
    const int idtemptemplocal = 9; // T*T
    const int idrhoh2local = 10; //rho*h2
    const int idrhoh2uzlocal = 11; //rho*h2*uz
    const int idrhoh2h2local = 12; //rho*h2*h2
    outNames[idrholocal] = "rho";
    outNames[iduzlocal] = "uz";
    outNames[idrhouzlocal] = "rho.uz";
    outNames[idrhouzuzlocal] = "rho.uz.uz";
    outNames[idrhohlocal] = "rho.h";
    outNames[idrhohuzlocal] = "rho.h.uz";
    outNames[idrhohhlocal] = "rho.h.h";
    outNames[idtemplocal] = "T";
    outNames[idtempuzlocal] = "T.uz";
    outNames[idtemptemplocal] = "T.T";
    outNames[idrhoh2local] = "rho.h2";
    outNames[idrhoh2uzlocal] = "rho.h2.uz";
    outNames[idrhoh2h2local] = "rho.h2.h2";
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<MultiFab*> indata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    const int nGrow = 0;
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const Vector<Real>& delta = amrData.DxLevel()[lev];
      const DistributionMapping dm(ba);
      indata[lev] = new MultiFab(ba,dm,nCompIn,nGrow);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(*indata[lev],lev,varNames,destFillComps); //Problem
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,0,&(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
        Array4<Real> const& outbox = (*outdata[lev]).array(mfi);
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  outbox(i,j,k,idrholocal) = inbox(i,j,k,idrholocal);
	  outbox(i,j,k,iduzlocal) = inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhouzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhouzuzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,iduzlocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhohlocal) = inbox(i,j,k,idrhohlocal);
	  outbox(i,j,k,idrhohuzlocal) = inbox(i,j,k,idrhohlocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhohhlocal) = inbox(i,j,k,idrhohlocal)*inbox(i,j,k,idrhohlocal)/inbox(i,j,k,idrholocal);
	  outbox(i,j,k,idtemplocal) = inbox(i,j,k,idtemplocal);
	  outbox(i,j,k,idtempuzlocal) = inbox(i,j,k,idtemplocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idtemptemplocal) = inbox(i,j,k,idtemplocal)*inbox(i,j,k,idtemplocal);
	  outbox(i,j,k,idrhoh2local) = inbox(i,j,k,idrholocal)*inbox(i,j,k,idh2local);
	  outbox(i,j,k,idrhoh2uzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,idh2local)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhoh2h2local) = inbox(i,j,k,idrholocal)*inbox(i,j,k,idh2local)*inbox(i,j,k,idh2local);
        });
	
      }
      
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_conservedQ");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
