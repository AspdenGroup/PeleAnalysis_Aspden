#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AppendToPlotFile.H>

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
  std::string names[nVars];
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
  if (pp.countval("probLo")==3) {
    for (int i=0; i<3; i++) {
      pp.get("probLo",probLo[i],i);
    }
  } else {
    amrex::Error("Need 3 values for probLo");
  }
  if (verbose) std::cout << "probLo = "
			 << probLo[0] << " " << probLo[1] << " " << probLo[2] << std::endl;
  
  Real probHi[3];
  if (pp.countval("probHi")==3) {
    for (int i=0; i<3; i++) {
      pp.get("probHi",probHi[i],i);
    }
  } else {
    amrex::Error("Need 3 values for probHi");
  }
  if (verbose) std::cout << "probHi = "
			 << probHi[0] << " " << probHi[1] << " " << probHi[2] << std::endl;

  Real dx[3];
  for (int i=0; i<3; i++) {
    dx[i] = (probHi[i]-probLo[i])/(Real)nx[i];
  }
  if (verbose) std::cout << "dx = "
			 << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;

  //
  // Let's try a slab domain decomposition
  //
  
  IntVect     pdLo(0,0,0);
  IntVect     pdHi(nx[0]-1,nx[1]-1,nx[2]-1);
  Box         probDomain(pdLo,pdHi);

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

  // make desination multifab (broken)
  MultiFab *mf;
  mf = new MultiFab(domainBoxArray,nVars,0);

  //
  // Populate data
  //
  if (verbose)
    std::cout << "Populating data:" << std::endl;

  // how do we loop over the boxes in the multifab properly?
  int iFile=0;
  for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {

    // destination fab
    FArrayBox& myFab = (*mf)[mfi];
    
    // load data
    std::cout << "iFile = " << iFile << ": " << infile[iFile].c_str() << std::endl;
    ifs.open(infile[iFile].c_str());
    fab.readFrom(ifs);
    ifs.close();
    const Box& inBox = fab.box();
    for (int dir=0; dir<2; dir++) {
      if (inBox.length(dir)!=nx[dir]) {
	std::cerr << "file = " << infile[iFile] << std::endl;
	std::cerr << "inBox.length(" << dir << ") = " << inBox.length(dir) << std::endl;
	amrex::Error("inBox.length mismatch!");
      }
    }

    // copy data
    myFab.copy(fab);
    
    iFile++;
  }

  // write the output plotfile
  // should be able to replace with modern call to writeplotfile

  if (verbose) {
    std::cout << "*** writing plotfile " << std::endl;
  }
  if(ParallelDescriptor::IOProcessor()) {
    if( ! amrex::UtilCreateDirectory(outfile, 0755)) {
      amrex::CreateDirectoryFailed(outfile);
    }
  }
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();
	   
  std::string outfileHeader(outfile);
  outfileHeader += "/Header";
  
  VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
  
  std::ofstream os;
  
  os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  
  int i;
  if(ParallelDescriptor::IOProcessor()) {
    if (verbose) {
	std::cout << "Opening file = " << outfileHeader << '\n';
    }
    
    os.open(outfileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if(os.fail()) {
      amrex::FileOpenFailed(outfileHeader);
    }
    //
    // Start writing plotfile.
    //
    os << "NavierStokes-V1.1" << '\n';
    os << nVars << '\n';
    for (int iVar=0; iVar<nVars; iVar++)
      os << names[iVar] << '\n';
    os << BL_SPACEDIM << '\n';
    os << time << '\n';
    int nLevels(1);
    os << nLevels-1 << '\n';
    for(i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for(i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    os << ' '; // refRatio
    os << '\n';
    os << probDomain;
    os << '\n';
    os << "0\n"; // level steps
    for(int k = 0; k < BL_SPACEDIM; k++) os << dx[k] << ' ';
    os << '\n';
    os << "0\n"; // coordSys
    os << "0\n"; // The bndry data width.
  }  // end if(ioproc)
	   
  //
  // Write out level by level.
  //
  int iLevel=0;
  int nGrids(domainBoxArray.size());
  char buf[64];
  sprintf(buf, "Level_%d", iLevel);
    
  if(ParallelDescriptor::IOProcessor()) {
    os << iLevel << ' ' << nGrids << ' ' << time << '\n';
    os << 0 << '\n';
      
    for(i = 0; i < nGrids; ++i) {
      for(int n = 0; n < BL_SPACEDIM; n++) {
	os << probLo[n] 
	   << ' '
	   << probHi[n]
	   << '\n';
      }
    }
    //
    // Build the directory to hold the MultiFabs at this level.
    //
    std::string Level(outfile);
    Level += '/';
    Level += buf;
      
    if( ! amrex::UtilCreateDirectory(Level, 0755)) {
      amrex::CreateDirectoryFailed(Level);
    }
  }  // end if(ioproc)
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();
    
  // Now build the full relative pathname of the MultiFab.
  static const std::string MultiFabBaseName("/MultiFab");
    
  std::string PathName(outfile);
  PathName += '/';
  PathName += buf;
  PathName += MultiFabBaseName;

  if (verbose) std::cout << "Writing data..." << std::endl;
  
  if(ParallelDescriptor::IOProcessor()) {
    // The full name relative to the Header file.
    std::string RelativePathName(buf);
    RelativePathName += MultiFabBaseName;
    os << RelativePathName << '\n';
  }
  VisMF::Write(*mf, PathName);
  
  os.close();
  
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "*** done writing plotfile " << std::endl;
  }
  DataServices::Dispatch(DataServices::ExitRequest, NULL);
  // ^^^ this calls ParallelDescriptor::EndParallel() and exit()

}
















