#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>

using namespace amrex;

//using namespace analysis_util;

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
    int nVar = 2+AMREX_SPACEDIM;
    Vector<std::string> varNames(nVar);
    varNames[0] = "RhoRT";
    varNames[1] = "density";
    varNames[2] = "x_velocity";
    varNames[3] = "y_velocity";
#if AMREX_SPACEDIM==3
    varNames[4] = "z_velocity";
#endif
    Real gamma=1.4;
    pp.query("gamma",gamma);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Vector<int> idVin(nVar);
    for (int i = 0; i<nVar; ++i) {
      idVin[i] = -1;
    } 
    
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));
    int coord = 0;
    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    

    for (int i=0; i<plotVarNames.size(); ++i)
    {
      for (int j = 0; j<nVar; ++j) {	
	if (plotVarNames[i] == varNames[j]) idVin[j] = i;
      }
    }
    for (int i = 0; i<nVar; ++i) {
      if (idVin[i] < 0) {
	//Print() << "Cannot find " << varNames[i] << " in pltfile data" << std::endl;
	Abort("Cannot find "+varNames[i]+" in pltfile data"); 
      }
    }
    const int nCompIn  = nVar;
    const int nCompOut = 1;
    Vector<std::string> outNames(1);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    for (int i = 0; i<nVar; i++) {
      destFillComps[i] = i;
      inNames[i] = varNames[i];
    }
    outNames[0] = "machNumber";
    
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    const int nGrow = 1;
    
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].define(ba,dm,nCompOut,nGrow);
      MultiFab indata(ba,dm,nCompIn,nGrow);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //Problem
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& varsIn  = indata.array(mfi);
        Array4<Real> const& varsOut = outdata[lev].array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  Real maxVel = 0;
	  for (int d = 0; d<AMREX_SPACEDIM; d++) {
	    maxVel += varsIn(i,j,k,2+d)*varsIn(i,j,k,2+d);
	  }
	  maxVel = std::sqrt(maxVel);
	  varsOut(i,j,k,0) = maxVel/std::sqrt(gamma*varsIn(i,j,k,0)/varsIn(i,j,k,1));
	});
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_mach");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms,0.0,isteps,refRatios);
  }
  Finalize();
  return 0;
}
