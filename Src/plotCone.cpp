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
    std::string fuelName="H2"; pp.query("fuelName",fuelName);
    
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
    Real r0,grad;
    pp.get("r0",r0);
    pp.get("grad",grad);
    int nCompIn = pp.countval("vars");
    const int nCompOut = nCompIn+1;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    if (nCompIn > 0) {
      pp.getarr("vars",inNames);
    }
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
      outNames[i] = inNames[i];
    }
    outNames[nCompIn] = "cone";
    Vector<Real> plo =  amrData.ProbLo();
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;
    
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);
      Real dx = amrData.DxLevel()[lev][0];
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //Problem
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = indata.array(mfi);
        Array4<Real> const& outbox = (*outdata[lev]).array(mfi);
    
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  for (int n = 0; n < nCompIn; n++) {
	    outbox(i,j,k,n) = inbox(i,j,k,n);
	  }
	  Real x = plo[0] + dx*(i+0.5);
	  Real y = plo[1] + dx*(j+0.5);
	  Real r = std::sqrt(x*x+y*y);
	  Real z = plo[2] + dx*(k+0.5);
	  Real p = r0+grad*z;
	  Real alpha = (r-p)/dx;
	  if (alpha < -0.5) {
	    outbox(i,j,k,nCompIn) = 1;
	  } else if (alpha > 0.5) {
	    outbox(i,j,k,nCompIn) = 0;
	  } else {
	    outbox(i,j,k,nCompIn) = 0.5-alpha;
	  }
	  
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_cone");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
