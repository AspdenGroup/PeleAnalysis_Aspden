#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
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
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

void
findNearestPointOnCenterline(Vector<Real> &x, Vector<Real> alpha, Vector<Real> coeffs, int dir, int otherdir, int slicedir) {
  if (alpha[otherdir] == 0.0) {
    x[otherdir] = 0.0;
    x[slicedir] = 0.0;
    x[dir] = 0.0;
    return;
  }
  Real tol = 1e-8;
  Real xmax = 0.008;
  Real xmin = 0.0;
  Real initial_guess = (xmax+xmin)/2.0;
  int itermax = 100;
  Real xn, xnp1,fxn, fxmin, fxmax;
  
  xn = initial_guess;
  fxn = std::pow(xn,2-coeffs[1])-alpha[dir]*std::pow(xn,1-coeffs[1]) + coeffs[0]*coeffs[0]*coeffs[1]*std::pow(xn,coeffs[1]) - alpha[otherdir]*coeffs[0]*coeffs[1];
  //std::cout << "fx0= " << fxn << std::endl;
  fxmin = std::pow(xmin,2-coeffs[1])-alpha[dir]*std::pow(xmin,1-coeffs[1]) + coeffs[0]*coeffs[0]*coeffs[1]*std::pow(xmin,coeffs[1]) - alpha[otherdir]*coeffs[0]*coeffs[1];
  //std::cout << "fxmin0= " << fxmin << std::endl;
  fxmax = std::pow(xmax,2-coeffs[1])-alpha[dir]*std::pow(xmax,1-coeffs[1]) + coeffs[0]*coeffs[0]*coeffs[1]*std::pow(xmax,coeffs[1]) - alpha[otherdir]*coeffs[0]*coeffs[1];
  //std::cout << "fxmax0= " << fxmax << std::endl;
  if(fxmin*fxmax > 0) {
    Abort("same sign");
  }
  if(fxn*fxmin > 0) {
    xmin = xn;
  } else {
    xmax = xn;
  }
  //std::cout << "xmin0= " << xmin << std::endl;
  //std::cout << "xmax0= " << xmax << std::endl;
  xnp1 = (xmax+xmin)/2.0;
     
  for (int i = 2; i<itermax; i++) {
    if(std::abs(xmax-xmin) > tol) {      
      xn = xnp1;
      fxn = std::pow(xn,2-coeffs[1])-alpha[dir]*std::pow(xn,1-coeffs[1]) + coeffs[0]*coeffs[0]*coeffs[1]*std::pow(xn,coeffs[1]) - alpha[otherdir]*coeffs[0]*coeffs[1];  
      //std::cout << "fx" << i-1 << "= " << fxn << std::endl;
      fxmin = std::pow(xmin,2-coeffs[1])-alpha[dir]*std::pow(xmin,1-coeffs[1]) + coeffs[0]*coeffs[0]*coeffs[1]*std::pow(xmin,coeffs[1]) - alpha[otherdir]*coeffs[0]*coeffs[1];
      //std::cout << "fxmin" << i-1 << "= " << fxmin << std::endl;
      if(fxn*fxmin > 0) {
	xmin = xn;
      } else {
	xmax = xn;
      }      
      //std::cout << "xmin" << i-1 << "= " << xmin << std::endl;
      //std::cout << "xmax" << i-1 << "= " << xmax << std::endl;
      xnp1 = (xmax+xmin)/2.0;      
      //std::cout << "x" << i << "= " << xnp1 << std::endl;
      
    } else {
      break;
    }
  }
  x[slicedir] = 0;
  x[dir] = xnp1;
  x[otherdir] = coeffs[0]*std::pow(xnp1,coeffs[1]);
  return;

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
    std::string maxComp="Y(H2)"; pp.query("maxComp",maxComp);
    int nVars;
    pp.countval("vars", nVars);
    if (nVars > 0) {
      Vector<std::string> vars; 
      vars.resize(nVars);
      pp.getarr("vars",vars);
    }
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

    int idMaxin = -1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == maxComp) idMaxin = i;
    }
    if (idMaxin<0)
      Print() << "Cannot find required data in pltfile" << std::endl;
    const int nCompIn  = 1;
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    inNames[0] = maxComp;
    destFillComps[0] = 0;
    
    int dir = 2; pp.query("dir",dir);
    int slicedir = 1; pp.query("slicedir",slicedir);    
    Box probDomain = amrData.ProbDomain()[finestLevel];
    Vector<Real> probLo = amrData.ProbLo();
    Vector<Real> probHi = amrData.ProbHi();
    Vector<Real> dx = amrData.DxLevel()[finestLevel];
    int numBoxes = probDomain.length(dir);
    int slice = (int)((numBoxes/2)-1); pp.query("slice",slice);    
    Box boxVector[numBoxes];
    AMREX_ASSERT(dir!=slicedir);
    bool foundOtherDir = false;
    int otherdir=0;
    while(!foundOtherDir) {
      if (otherdir!=dir && otherdir!=slicedir) {
	foundOtherDir = true;
      } else {
	otherdir += 1;
      }      
    }
    for (int i = 0; i < numBoxes; i++) {      
      IntVect boxBottom;
      IntVect boxTop;
      boxBottom[dir] = i;
      boxBottom[slicedir] = slice;
      boxBottom[otherdir] = 0;
      boxTop[dir] = i;
      boxTop[slicedir] = slice;
      boxTop[otherdir] = numBoxes-1;
      boxVector[i] = Box(boxBottom,boxTop);
    }
    const BoxArray centerlineBA = BoxArray(boxVector,numBoxes);
    
    MultiFab CLMF(centerlineBA,DistributionMapping(centerlineBA),nCompIn,0);
    Vector<Real> x,y,H2val;
    x.resize(numBoxes);
    y.resize(numBoxes);
    H2val.resize(numBoxes);
    Print() << "Reading data" << std::endl;
    amrData.FillVar(CLMF,finestLevel,inNames,destFillComps); //Problem
    Print() << "Data has been read" << std::endl;  
    for (MFIter mfi(CLMF,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& varsIn = CLMF.array(mfi);
      Real compVal = 0.0;
      IntVect maxLoc = {0,0,0};
      AMREX_PARALLEL_FOR_3D (bx, i,  j, k, {
	  if (varsIn(i,j,k,0) > compVal) {
	    maxLoc[0] = i;
	    maxLoc[1] = j;
	    maxLoc[2] = k;
	    compVal = varsIn(i,j,k,0);
	  }
	});
      x[mfi.index()] = probLo[dir] + dx[dir]*(0.5+maxLoc[dir]);
      y[mfi.index()] = probLo[otherdir] + dx[otherdir]*(0.5+maxLoc[otherdir]);
      H2val[mfi.index()] = compVal;
    }
    ParallelDescriptor::ReduceRealSum(x.dataPtr(),numBoxes,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(y.dataPtr(),numBoxes,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(H2val.dataPtr(),numBoxes,ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
      std::string outfile = "centerline.dat";
      pp.query("outfile",outfile);
      // Write ASCII output file
      std::ofstream os(outfile.c_str(),std::ios::out);
      for (int n = 0; n < numBoxes; n++) {
	os << x[n] << " " << y[n] << " " << H2val[n] << std::endl;
      }
      os.close();
    }
    int nCoeffs = pp.countval("coeffs");
    if (nCoeffs > 0) {
      AMREX_ASSERT(nCoeffs==2);
      Vector<Real> coeffs;
      coeffs.resize(nCoeffs);
      pp.getarr("coeffs",coeffs);
      const int nCompOut = 3;
      Vector<std::string> outNames(nCompOut);
      outNames[0] = "s";
      outNames[1] = "r";
      outNames[2] = "theta";
      const BoxArray originalba = amrData.boxArray(finestLevel);
      
      MultiFab outdata(originalba,DistributionMapping(originalba),nCompOut,0);
    
      for (MFIter mfi(outdata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  Array4<Real> const& array = outdata.array(mfi);
	  
	  AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
				{
				  Vector<Real> alpha;
				  alpha.resize(3);
				  alpha[0] = probLo[0] + dx[0]*(0.5+i);
				  alpha[1] = probLo[1] + dx[1]*(0.5+j);
				  alpha[2] = probLo[2] + dx[2]*(0.5+k);
				  Vector<Real> x;
				  x.resize(3);
				  findNearestPointOnCenterline(x,alpha,coeffs,dir,otherdir,slicedir);
				  //Print() << "(" << alpha[0] << ", " << alpha[1] << ", " << alpha[2] << ") -> (" << x[0] << ", " << x[1] << ", " << x[2]  << ")" << std::endl;
				  array(i,j,k,0) = std::sqrt(x[dir]*x[dir] + x[otherdir]*x[otherdir]);
				  Vector<Real> xalpha;
				  xalpha.resize(3);
				  xalpha[0] = alpha[0]-x[0];
				  xalpha[1] = alpha[1]-x[1];
				  xalpha[2] = alpha[2]-x[2];
				  array(i,j,k,1) = std::sqrt(xalpha[0]*xalpha[0] + xalpha[1]*xalpha[1] + xalpha[2]*xalpha[2]);
				  array(i,j,k,2) = 0;
				  if (array(i,j,k,1) > 0)
				    array(i,j,k,2) = std::acos(xalpha[slicedir]/(array(i,j,k,1)));
				  
				});
	  //std::cout << mfi.index() << std::endl;
	}
      RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
      Vector<int> is_per(AMREX_SPACEDIM,0);
      Geometry geoms(amrData.ProbDomain()[finestLevel],&rb,0,&(is_per[0]));
      
      std::string outfile(getFileRoot(plotFileName) + "_coords");
      Print() << "Writing new data to " << outfile << std::endl;
      const bool verb = false;
      int level_step;
      WriteSingleLevelPlotfile(outfile,outdata,outNames,geoms,0,level_step);
      
    }
      
  }
  Finalize();
  return 0;
}

  
