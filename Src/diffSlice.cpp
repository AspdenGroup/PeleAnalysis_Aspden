#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);

  int nProcs = ParallelDescriptor::NProcs();
  if (nProcs>1)
    amrex::Error("Not yet fully-implemented in parallel!");

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
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  // read in list of fabs to load
  int nFiles=pp.countval("infile");
  std::string infile[nFiles];
  if (verbose) std::cout << "infiles =";
  for (int iFile=0; iFile<nFiles; iFile++) {
    pp.get("infile",infile[iFile],iFile);
    if (verbose) std::cout << " " << infile[iFile]; 
  }
  if (verbose) std::cout << std::endl;

  // read plot file name to output
  std::string outfile;
  pp.get("outfile",outfile);
  if (verbose) std::cout << "outfile = " << outfile << std::endl;

  // read a time to assign
  Real time(0.);
  pp.query("time",time);
  if (verbose) std::cout << "time = " << time << std::endl;

  // read in the variable names to use
  int nVars=pp.countval("names");
  Vector<std::string> names; names.resize(nVars);
  if (verbose) std::cout << "variables =";
  for (int iVar=0; iVar<nVars; iVar++) {
    pp.get("names",names[iVar],iVar);
    if (verbose) std::cout << " " << names[iVar]; 
  }
  if (verbose) std::cout << std::endl;

  //
  // Read fab data (only using header at this stage)
  //
  if (verbose) std::cout << "Reading fab..." << std::endl;
  std::ifstream ifs;
  ifs.open(infile[0].c_str());
  FArrayBox fab;
  fab.readFrom(ifs);
  ifs.close();
  if (verbose) std::cout << "   ... done." << std::endl;
  
  std::cout << "fab.nComp = " << fab.nComp() << std::endl;
  if (fab.nComp()!=nVars) {
    amrex::Error("Mismatch fab.nComp != nVars");
  }

  // figure out how big the full domain box needs to be
  const Box& box = fab.box();
  std::cout << "fab.box().length = " << box.length() << std::endl;

  // number of cells in each direction
  // (assume stacking in z)
  int nx[3];
  nx[0] = box.length(0);
  nx[1] = box.length(1);
  nx[2] = nFiles;
  int pCells=nx[0]*nx[1];
  int nCells=nx[0]*nx[1]*nx[2];
    
  if (verbose) std::cout << "cells = "
			 << nx[0] << " " << nx[1] << " " << nx[2] << " "
			 << " (" << nCells << ")" << std::endl;

  // assign a physical size
  Real probLo[3];
  int calcz;
  if (pp.countval("probLo")==3) {
    for (int i=0; i<3; i++) {
      pp.get("probLo",probLo[i],i);
    }
    calcz = 0; 
  } else if (pp.countval("probLo") == 2) {
    for (int i=0; i<2; i++) {
      pp.get("probLo",probLo[i],i);
    }
    calcz = 1;
  } else {
    amrex::Error("Need 3 values for probLo");
  }

  
  Real probHi[3];
  if (pp.countval("probHi")==3) {
    for (int i=0; i<3; i++) {
      pp.get("probHi",probHi[i],i);
    }
  } else {
    amrex::Error("Need 3 values for probHi");
  }
  Vector<int> is_per(AMREX_SPACEDIM,0); //hard code to no periodicity for the minute
  Real dx[3];
  if (calcz) {
    for (int i=0; i<2; i++) {
      dx[i] = (probHi[i]-probLo[i])/(Real)nx[i];
    }
    dx[2] = dx[0];
    probLo[2] = -dx[2]*nx[2];
  } else {
    for (int i=0; i<3; i++) {
      dx[i] = (probHi[i]-probLo[i])/(Real)nx[i];
    }
  }

  if (verbose) std::cout << "probLo = "
			 << probLo[0] << " " << probLo[1] << " " << probLo[2] << std::endl;
  
  if (verbose) std::cout << "probHi = "
			 << probHi[0] << " " << probHi[1] << " " << probHi[2] << std::endl;
  
  if (verbose) std::cout << "dx = "
			 << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;
  RealBox rb(probLo,probHi); // make real box for geometry
  
  //
  // Let's try a slab domain decomposition
  //
  
  IntVect     pdLo(0,0,0);
  IntVect     pdHi(nx[0]-1,nx[1]-1,nx[2]-1);
  Box         probDomain(pdLo,pdHi);
  int coord = 0; //hard code cartesian
  Geometry geoms(probDomain, &rb, coord, &(is_per[0]));
  // make a box for each file
  int         nBoxes(nFiles);
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
  
  // And now the distibution mapping
  Vector<int> pmap(nBoxes);
  for (int iBox=0; iBox<nBoxes; iBox++)
    pmap[iBox] = myProc; // maybe think about hacking this to parallelise this step
  DistributionMapping domainDistMap(pmap);

  
  MultiFab *mf;
  mf = new MultiFab(domainBoxArray,domainDistMap,nVars,0);
  //
  // Populate data
  //
  if (verbose)
    std::cout << "Populating data:" << std::endl;

  // how do we loop over the boxes in the multifab properly?
  int iFile=0;
  IntVect topslicebig = pdHi;
  IntVect topslicesmall = pdLo;
  topslicebig[2] = pdLo[2];
  Box topBox(topslicesmall,topslicebig);
  FArrayBox topFab(topBox,nVars);
  int slicenum = 0;

  ifs.open(infile[slicenum].c_str());
  fab.readFrom(ifs);
  ifs.close();
  IntVect shift = {0,0,-fab.box().smallEnd(2)};
  fab.shift(shift);  
  topFab.copy(fab);
    
  int diffComp=3;
  int corrComp=diffComp;
  int pnorm=1;
  FArrayBox topFabsq(topBox,nVars);
  topFabsq.copy(topFab);
  topFabsq.minus(topFab.sum(corrComp)/(Real)pCells);
  topFabsq.mult(topFabsq);
 
  FArrayBox diffFab(topBox,nVars);
  Vector<Real> diffs(nFiles);
  FArrayBox corrFab(topBox,nVars);
  Vector<Real> autocorr(nFiles);
  for (iFile=0; iFile<nFiles; iFile++) {
    
    // load data
    std::cout << "iFile = " << iFile << ": " << infile[iFile].c_str() << std::endl;
    ifs.open(infile[iFile].c_str());
    fab.readFrom(ifs);
    ifs.close();
    
    IntVect shift = {0,0,-fab.box().smallEnd(2)};
    fab.shift(shift);
    
    diffFab.copy(topFab);
    diffFab.minus(fab);
    diffs[iFile] = diffFab.norm(pnorm,diffComp,1);

    corrFab.copy(topFab);
    corrFab.minus(corrFab.sum(corrComp)/(Real)pCells);
    fab.minus(fab.sum(corrComp)/(Real)pCells);
    corrFab.mult(fab);
    autocorr[iFile] = corrFab.sum(corrComp)/topFabsq.sum(corrComp);
    //    iFile++;
  }

  

  // write the output plotfile
  // should be able to replace with modern call to writeplotfile

  if (verbose) {
    std::cout << "*** writing plotfile " << std::endl;
  }
  //  int levelSteps;
  //WriteSingleLevelPlotfile(outfile,*mf,names, geoms,time,levelSteps);

  outfile="diff.dat";
  std::ofstream os(outfile.c_str(),std::ios::out);
  for (int n = 0; n<nFiles; n++) {
    os << n << " " << diffs[n] << " " << autocorr[n] << "\n";
  }
  os.close();
}
















