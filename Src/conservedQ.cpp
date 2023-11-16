#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

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
    int avgDown = 0; pp.query("avgDown",avgDown);
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
    Vector<std::string> varNames = {"density","z_velocity","rhoh","temp","avg_pressure","Y(H2)","x_velocity","y_velocity"};

    const int nCompIn  = varNames.size();
    int nCompOut = 21;
    
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
    const int idpressurelocal = 4; //avg_pressure - keep
    const int idh2in = 5; //h2 (just for rhoh2)
    const int idxvelin = 6; //xvel (just for ur)
    const int idyvelin = 7; //yvel (just for ur)
    //out
    const int idrhouzlocal = 5; // rho*uz
    const int idrhouzuzlocal = 6; // rho*uz*uz
    const int idurlocal = 7; // ur
    const int idrhourlocal = 8; // rho*ur
    const int idrhoururlocal = 9; // rho*ur*ur
    const int idrhouzurlocal = 10; // rho*ur*uz
    const int idrhohuzlocal = 11; //rho*h*uz
    const int idrhohhlocal = 12; //rho*h*h
    const int idrhohurlocal = 13; //rho*h*ur
    const int idtempuzlocal = 14; // T*uz
    const int idtempurlocal = 15; //T*ur
    const int idtemptemplocal = 16; // T*T
    const int idrhoh2local = 17; //rho*h2
    const int idrhoh2uzlocal = 18; //rho*h2*uz
    const int idrhoh2urlocal = 19; //rho*h2*ur
    const int idrhoh2h2local = 20; //rho*h2*h2
    outNames[idrholocal] = "rho";
    outNames[iduzlocal] = "uz";
    outNames[idrhouzlocal] = "rho.uz";
    outNames[idrhouzuzlocal] = "rho.uz.uz";
    outNames[idrhouzurlocal] = "rho.uz.ur";
    outNames[idurlocal] = "ur";
    outNames[idrhourlocal] = "rho.ur";
    outNames[idrhoururlocal] = "rho.ur.ur";
    outNames[idrhohlocal] = "rho.h";
    outNames[idrhohuzlocal] = "rho.h.uz";
    outNames[idrhohurlocal] = "rho.h.ur";
    outNames[idrhohhlocal] = "rho.h.h";
    outNames[idtemplocal] = "T";
    outNames[idtempuzlocal] = "T.uz";
    outNames[idtempurlocal] = "T.ur";
    outNames[idtemptemplocal] = "T.T";
    outNames[idrhoh2local] = "rho.h2";
    outNames[idrhoh2uzlocal] = "rho.h2.uz";
    outNames[idrhoh2urlocal] = "rho.h2.ur";
    outNames[idrhoh2h2local] = "rho.h2.h2";
    outNames[idpressurelocal] = "pressure";
    Vector<MultiFab*> outdata(Nlev);
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
      outdata[lev] = new MultiFab(ba,dm,nCompOut,nGrow);
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(*indata[lev],lev,varNames,destFillComps); //Problem
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,0,&(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      Vector<Real> dx = amrData.DxLevel()[lev];
      Vector<Real> plo = amrData.ProbLo();
      
      for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
        Array4<Real> const& outbox = (*outdata[lev]).array(mfi);
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  Real x = plo[0] + (i+0.5)*dx[0];
	  Real y = plo[1] + (j+0.5)*dx[1];
	  Real r = std::sqrt(x*x+y*y);
	  outbox(i,j,k,idrholocal) = inbox(i,j,k,idrholocal);
	  outbox(i,j,k,iduzlocal) = inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idurlocal) = (x*inbox(i,j,k,idxvelin)+y*inbox(i,j,k,idyvelin))/r;
	  outbox(i,j,k,idrhouzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhouzuzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,iduzlocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhouzurlocal) = inbox(i,j,k,idrholocal)*outbox(i,j,k,idurlocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhourlocal) = inbox(i,j,k,idrholocal)*outbox(i,j,k,idurlocal);
	  outbox(i,j,k,idrhoururlocal) = inbox(i,j,k,idrholocal)*outbox(i,j,k,idurlocal)*outbox(i,j,k,idurlocal);
	  outbox(i,j,k,idrhohlocal) = inbox(i,j,k,idrhohlocal);
	  outbox(i,j,k,idrhohuzlocal) = inbox(i,j,k,idrhohlocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhohurlocal) = inbox(i,j,k,idrhohlocal)*outbox(i,j,k,idurlocal);
	  outbox(i,j,k,idrhohhlocal) = inbox(i,j,k,idrhohlocal)*inbox(i,j,k,idrhohlocal)/inbox(i,j,k,idrholocal);
	  outbox(i,j,k,idtemplocal) = inbox(i,j,k,idtemplocal);
	  outbox(i,j,k,idtempuzlocal) = inbox(i,j,k,idtemplocal)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idtempurlocal) = inbox(i,j,k,idtemplocal)*outbox(i,j,k,idurlocal);
	  outbox(i,j,k,idtemptemplocal) = inbox(i,j,k,idtemplocal)*inbox(i,j,k,idtemplocal);
	  outbox(i,j,k,idrhoh2local) = inbox(i,j,k,idrholocal)*inbox(i,j,k,idh2in);
	  outbox(i,j,k,idrhoh2uzlocal) = inbox(i,j,k,idrholocal)*inbox(i,j,k,idh2in)*inbox(i,j,k,iduzlocal);
	  outbox(i,j,k,idrhoh2urlocal) = outbox(i,j,k,idrhoh2local)*outbox(i,j,k,idurlocal);
	  outbox(i,j,k,idrhoh2h2local) = outbox(i,j,k,idrhoh2local)*inbox(i,j,k,idh2in);
	  outbox(i,j,k,idpressurelocal) = inbox(i,j,k,idpressurelocal);
        });
	
      }
      
      Print() << "Derive finished for level " << lev << std::endl;
    }
    std::string outfile(getFileRoot(plotFileName) + "_conservedQ");
    const bool verb = false;
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    if (avgDown) {
      Print() << "Averaging down to coarse and writing to "+outfile << std::endl;
      for (int lev = finestLevel; lev > 0; --lev) {
	MultiFab* fineMF = outdata[lev];
	MultiFab* coarseMF = outdata[lev-1];
	average_down(*fineMF,*coarseMF,0,nCompOut,refRatios[lev-1]);
      }
      amrex::WriteSingleLevelPlotfile(outfile,*outdata[0],outNames,geoms[0],0.0,0);
    } else {
      Print() << "Writing to " +outfile << std::endl;
      Vector<int> isteps(Nlev, 0);
      amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms, 0.0, isteps, refRatios);
    }

  }
  Finalize();
  return 0;
}
