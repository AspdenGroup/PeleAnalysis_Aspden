#include <string>
#include <iostream>

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
  std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
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

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Set up input field data names, and destination components to load data upon read.
    int ncomps = pp.countval("keepVars");
    Vector<std::string> inNames(ncomps);
    Vector<int> destFillComps(ncomps);
    Print() << "Keeping vars: " << std::endl;
    for (int i = 0; i<ncomps; i++) {
      pp.get("keepVars",inNames[i],i);
      destFillComps[i] = i;
      Print() << inNames[i] << std::endl;
    }
    
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();

    Vector<int> idVars(ncomps);
    for (int i=0; i<ncomps; i++) {
      idVars[i] = -1;
      for (int j=0; j<plotVarNames.size(); j++) {
	if (plotVarNames[j] == inNames[i]) {
	  idVars[i] = j;
	  break;
	}
      }
    }
    for (int i=0; i<ncomps; i++) {
      if (idVars[i] < 0) {
	Print() << inNames[i] << " not found" << std::endl;
      }
    }
    
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    Vector<std::string> outNames = inNames;

    // Loop over AMR levels in the plotfile, read the data and do work
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    const int nGrow = 0;
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    
    int coord = 0;
    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      outdata[lev] =  MultiFab(ba,dm,outNames.size(),nGrow); 
      Print() << "Loading in level " << lev << std::endl;
      amrData.FillVar(outdata[lev],lev,inNames,destFillComps);
    }
      
    std::string outfile(getFileRoot(plotFileName) + "_trim");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});

    WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);
    
  }
  Finalize();
  return 0;
}
