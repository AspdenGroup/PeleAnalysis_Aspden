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
    Vector<std::string> dummyVars;
    Vector<Real> dummyVals;
    int nDummyVars = pp.countval("dummyVars");
    AMREX_ALWAYS_ASSERT(pp.countval("dummyVals") == nDummyVars);
    if (nDummyVars > 0) {
      dummyVars.resize(nDummyVars);
      pp.getarr("dummyVars",dummyVars);
      dummyVals.resize(nDummyVars);
      pp.getarr("dummyVals",dummyVals);
    }
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const int nCompIn  = plotVarNames.size();
    const int nCompOut = nCompIn+nDummyVars;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
      inNames[i] = plotVarNames[i];
      outNames[i] = plotVarNames[i];
    }
    for (int i = 0; i<nDummyVars; i++) {
      outNames[nCompIn+i] = dummyVars[i];
    }
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;
    int b[3] = {1, 1, 1};
   

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
	Array4<Real> const& inarr  = indata.array(mfi);
        Array4<Real> const& outarr = (*outdata[lev]).array(mfi);
	
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  for (int n = 0; n<nCompIn; n++) {
	    outarr(i,j,k,n) = inarr(i,j,k,n);
	  }
	  for (int n= 0; n<nDummyVars;n++) {
	    outarr(i,j,k,n+nCompIn) = dummyVals[n];
	  }
	});
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_dummy");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
