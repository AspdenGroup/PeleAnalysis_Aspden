#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

void
writeSingleArray(Vector<Real> arr, int arrsize, std::string filename) {
  FILE *file;
  file = fopen(filename.c_str(),"w");
  for (int i = 0; i < arrsize; i++) {
    fprintf(file,"%e ",arr[i]);
  }
  fclose(file);
  return;
}


void
writeDoubleArray(Vector<Vector<Real>> arr, int arrsize1, int arrsize2, std::string filename) {
  FILE *file;
  file = fopen(filename.c_str(),"w");
  for (int i = 0; i < arrsize1; i++) {
    for (int j = 0; j < arrsize2; j++) {
      fprintf(file,"%e ",arr[i][j]);
    }
    fprintf(file,"\n");
  }
  fclose(file);
  return;
}


int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);

  int nProcs = ParallelDescriptor::NProcs();
  int myProc = ParallelDescriptor::MyProc();

  ParmParse pp;
  
  bool verbose(false);
  if(pp.contains("verbose") || (pp.contains("v"))) {
    verbose = true;
    AmrData::SetVerbose(true);
  }

  if (ParallelDescriptor::IOProcessor()) {
    verbose = true;
    AmrData::SetVerbose(true);
  }
  
  DataServices::SetBatchMode();
  // Open plotfile header and create an amrData object pointing into it
  std::string plotFileName; pp.get("infile",plotFileName);
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(plotFileName, fileType);
  if( ! dataServices.AmrDataOk()) {
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  AmrData& amrData = dataServices.AmrDataRef();

  // Set up input field data names, and destination components to load data upon read.
  int ncomps = pp.countval("vars");
  AMREX_ALWAYS_ASSERT(ncomps==2);
  Vector<std::string> vars(ncomps);
  Vector<int> destfillcomps(ncomps);
  Print() << "JPDF of : " << std::endl;
  for (int i = 0; i<ncomps; i++) {
    pp.get("vars",vars[i],i);
    destfillcomps[i] = i;
    Print() << vars[i] << std::endl;
  }
  
  const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
  
  Vector<int> idVars(ncomps);
  for (int i=0; i<ncomps; i++) {
    idVars[i] = -1;
    for (int j=0; j<plotVarNames.size(); j++) {
      if (plotVarNames[j] == vars[i]) {
	idVars[i] = j;
	break;
      }
    }
  }
  for (int i=0; i<ncomps; i++) {
    if (idVars[i] < 0) {
      Print() << vars[i] << " not found" << std::endl;
    }
  }

  RealBox rb(&(amrData.ProbLo()[0]), 
	     &(amrData.ProbHi()[0]));
  
  // Loop over AMR levels in the plotfile, read the data and do work
  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel + 1;
  const int nGrow = 0;
  int dir = 2;
  pp.query("dir",dir);
  Vector<MultiFab> indata(Nlev);
  
  //
  // Let's try a slab domain decomposition
  //

  for (int lev = 0; lev<Nlev; lev++) {  
  
  
    Box probDomain = amrData.ProbDomain()[lev];
    int nslices = probDomain.bigEnd(dir);
  
    // make a box for each file
    int         nBoxes(nslices);
    BoxArray    domainBoxArray(nBoxes);
    
    // Make some the slabs
    Box         tempBox(probDomain);    
    Vector<int> tempBoxSmall(nBoxes,0);
    Vector<int> tempBoxBig(nBoxes,0);
    
    if (verbose)
      std::cout << "Slab decomposition:" << std::endl;
    for (int iBox=0; iBox<nBoxes; iBox++) {
      tempBoxSmall[iBox] = probDomain.smallEnd(2) + iBox;
      tempBoxBig[iBox]   = tempBoxSmall[iBox];
      if (verbose) {
	std::cout << "   iBox / small / big: "
		  << iBox << " / "
		  << tempBoxSmall[iBox] << " / "
		  << tempBoxBig[iBox] << std::endl;
      }
    }    
    for (int iBox=0; iBox<nBoxes; iBox++) {
      tempBox.setSmall(Amrvis::ZDIR, tempBoxSmall[iBox]);
      tempBox.setBig(Amrvis::ZDIR, tempBoxBig[iBox]); 
      domainBoxArray.set(iBox, tempBox);
    }
    
    Vector<int> pmap(nBoxes);
    for (int iBox=0; iBox<nBoxes; iBox++)
      pmap[iBox] = myProc+iBox%nProcs; // maybe think about hacking this to parallelise this step
    DistributionMapping domainDistMap(pmap);

    indata[lev] = MultiFab(domainBoxArray,domainDistMap,ncomps,0);
    amrData.FillVar(indata[lev],lev,vars,destfillcomps);
  }
  Vector<Vector<Real>> dxArr = amrData.DxLevel();
  //
  int nBins = 128;
  Real minv1 = indata[0].min(0);
  Real minv2 = indata[0].min(1);
  Real maxv1 = indata[0].max(0);
  Real maxv2 = indata[0].max(1);
  pp.query("bins",nBins);
  pp.query("minv1",minv1);
  pp.query("minv2",minv2);
  pp.query("maxv1",maxv1);
  pp.query("maxv2",maxv2);
  Real dv1 = (maxv1-minv1)/(Real)nBins;
  Real dv2 = (maxv2-minv2)/(Real)nBins;
  Vector<Real> v1(nBins);
  Vector<Real> v2(nBins);
  for (int n = 0; n < nBins; n++) {
    v1[n] = minv1 + (n+0.5)*dv1;
    v2[n] = minv2 + (n+0.5)*dv2;
  }
  Vector<Real> v1tmp(nBins,0);
  Vector<Vector<Real>> v2tmp(nBins,v1tmp);
  
  Box probDomain0 = amrData.ProbDomain()[0];
  int nBoxes0 = probDomain0.bigEnd(dir);
  int ncells = probDomain0.bigEnd(0)*probDomain0.bigEnd(1);
  Vector<Vector<Vector<Real>>> JPDFs(nBoxes0,v2tmp);
  //just loop base grid for now - higher levels will require more thought
  Real dxCoarse = dxArr[0][0];
  int iSlice = 0;
  //also scatter plot
  Vector<Vector<Vector<Real>>> ZT(nBoxes0);
  for (MFIter mfi(indata[0]); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    Print() << bx << std::endl;
    Array4<Real> const& inarr = indata[0].array(mfi);
    Vector<Vector<Real>> ZTloc;
  
    AMREX_PARALLEL_FOR_3D (bx, i, j, k,
			   {
			     int v1idx = (int)((Real)nBins*(inarr(i,j,k,0)-minv1)/(maxv1-minv1));
			     int v2idx = (int)((Real)nBins*(inarr(i,j,k,1)-minv2)/(maxv2-minv2));
			     if (v1idx >= 0 && v2idx >= 0 && v1idx < nBins && v2idx < nBins) {
			       JPDFs[iSlice][v2idx][v1idx] += dxCoarse*dxCoarse;
			     }
			     ZTloc.push_back({inarr(i,j,k,0),inarr(i,j,k,1)});
			   });
    
    ZT[iSlice] = ZTloc;
    
    iSlice++;
  }
  
  std::string v1filename = plotFileName+"/"+vars[0]+".dat";
  std::string v2filename = plotFileName+"/"+vars[1]+".dat";
  writeSingleArray(v1,nBins,v1filename);
  writeSingleArray(v2,nBins,v2filename);
  for (int n=0; n<nBoxes0; n++) {
    std::string filename = plotFileName+"/"+"slice"+std::to_string(n)+"JPDF.dat";
    writeDoubleArray(JPDFs[n],nBins,nBins,filename);
    filename = plotFileName+"/"+"slice"+std::to_string(n)+"scatter.dat";
    writeDoubleArray(ZT[n],ncells,2,filename);
  }  
}

















