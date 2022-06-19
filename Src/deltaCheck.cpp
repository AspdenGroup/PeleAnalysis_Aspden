#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>
/*
#include <mechanism.h>
#include <chemistry_file.H>
#include <util.H>
#include <util_F.H>
#include <Transport_F.H>
#include <Fuego_EOS.H>
*/
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
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
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
    Vector<std::string> varNames;
    int nVar= pp.countval("vars");
    varNames.resize(nVar);
    pp.getarr("vars",varNames,0,nVar);
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
    const int nCompOut = nVar;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    for (int i = 0; i<nVar; i++) {
      destFillComps[i] = i;
      inNames[i] = varNames[i];
      outNames[i] = "maxDelta_"+varNames[i];
    }
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 1;
    //int b[3] = {1, 1, 1};
    
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //Problem
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& varsIn  = indata.array(mfi);
        Array4<Real> const& varsOut = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  for (int n=0;n<nVar; n++) {
	    Real deltaxp1 = std::abs(varsIn(i+1,j,k,n)-varsIn(i,j,k,n));
	    Real deltaxm1 = std::abs(varsIn(i,j,k,n)-varsIn(i-1,j,k,n));
	    Real deltayp1 = std::abs(varsIn(i,j+1,k,n)-varsIn(i,j,k,n));
	    Real deltaym1 = std::abs(varsIn(i,j,k,n)-varsIn(i,j-1,k,n));
	    Real deltazp1 = std::abs(varsIn(i,j,k+1,n)-varsIn(i,j,k,n));
	    Real deltazm1 = std::abs(varsIn(i,j,k,n)-varsIn(i,j,k-1,n));
	    varsOut(i,j,k,n) = std::max(std::max(deltaxp1,std::max(deltaxm1,deltayp1)),std::max(deltaym1,std::max(deltazp1,deltazm1)));
	  } 
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_delta");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = true;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
