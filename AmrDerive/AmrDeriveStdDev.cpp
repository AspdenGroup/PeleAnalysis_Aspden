
//
// This is a version of AmrDerive.cpp that calculates integrals of
// quantities and writes out scalars instead of plotfiles.
//

#include <new>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <unistd.h>

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "VisMF.H"
#include "Geometry.H"

#include "AmrDeriveStdDev_F.H"

// This one is for the ghost cells (filcc)
#include "xtra_F.H"


int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  if(argc == 1)
    PrintUsage(argv[0]);

  ParmParse pp;

  if(pp.contains("help"))
    PrintUsage(argv[0]);

  FArrayBox::setFormat(FABio::FAB_IEEE_32);

  bool verbose;
  if(pp.contains("verbose") || (pp.contains("v"))) {
    verbose = true;
    AmrData::SetVerbose(true);
  }

  // Let's set verbose to true, without setting AmrData(verbose)
  if (ParallelDescriptor::IOProcessor())
    verbose = true;
  else
    verbose = false;

  // Let's see how many processors we're working with here
  int nProcs(ParallelDescriptor::NProcs());    //VSHOWVAL(verbose, nProcs);
  int myProc = ParallelDescriptor::MyProc();

  // Count number of input plot files
  int nPlotFiles(pp.countval("infile"));
  if(nPlotFiles <= 0) {
    std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
    std::cerr << "Exiting." << std::endl;
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  //VSHOWVAL(verbose, nPlotFiles);

  // Make an array of srings containing paths of input plot files
  Array<std::string> plotFileNames(nPlotFiles);
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    pp.get("infile", plotFileNames[iPlot], iPlot);
    //VSHOWVAL(verbose, plotFileNames[iPlot]);
  }

  // More random initialisation stuff
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
  Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
  Array<Real>           time(nPlotFiles);

  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    if (verbose && ParallelDescriptor::IOProcessor())
      //std::cout << "Loading " << plotFileNames[iPlot] << std::endl;

    dataServicesPtrArray[iPlot] = new DataServices(plotFileNames[iPlot], fileType);               // Populate DataServices array

    if( ! dataServicesPtrArray[iPlot]->AmrDataOk())                                               // Check AmrData ok
      DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not

    amrDataPtrArray[iPlot] = &(dataServicesPtrArray[iPlot]->AmrDataRef());                        // Populate DataPtrArray

    time[iPlot] = amrDataPtrArray[iPlot]->Time();

    //if (verbose) //std::cout << "Time = " << time[iPlot] << std::endl;
  }

  int finestLevel = amrDataPtrArray[0]->FinestLevel();
  pp.query("finestLevel",finestLevel);
  int nLevels = finestLevel + 1;

  //Pick axis to integrate over
  int axis = 2;
  pp.query("axis",axis);
  //
  // Variables to load
  // may want to hard wire this bit
  //
  int nVars(pp.countval("vars"));
  Array<string> varName(nVars);
  for(int v = 0; v < nVars; v++) {
    pp.get("vars", varName[v], v);
  }
  Array<int> destFillComps(nVars);
  for (int i=0; i<nVars; ++i)
    destFillComps[i] = i;

  // only necessary to set to 1 if derivatives required
  int nGrow = 0;

  //
  // Loop over plots
  //
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

      //
      // Load the data
      //
      PArray<MultiFab> mf(nLevels,PArrayManage);
      for (int iLevel=0; iLevel<nLevels; iLevel++) {
	  const BoxArray& ba = amrDataPtrArray[iPlot]->boxArray(iLevel);
	  mf.set(iLevel,new MultiFab(ba,nVars+1,nGrow)); // Do nVars+1 so that we can use the last one as a flag
	  if (ParallelDescriptor::IOProcessor())
	      //std::cout << "Reading level: " << iLevel << std::endl;
	  amrDataPtrArray[iPlot]->FillVar(mf[iLevel],iLevel,varName,destFillComps);   // FillVar
	  if (ParallelDescriptor::IOProcessor())
	      //std::cout << "...done reading level: " << iLevel << std::endl;
	  for (int iVar=0; iVar<nVars; iVar++)
	      amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varName[iVar]));   // Flush grids
      }

      const Real *dx = amrDataPtrArray[iPlot]->DxLevel()[finestLevel].dataPtr();
      Box probDomain = amrDataPtrArray[iPlot]->ProbDomain()[finestLevel];

      int isize(probDomain.length(Amrvis::XDIR));
      int jsize(probDomain.length(Amrvis::YDIR));
      int ksize(probDomain.length(Amrvis::ZDIR));
      int size;
      if (axis == 0) {
        size = isize;
      } else if (axis == 1) {
        size = jsize;
      } else if (axis == 2) {
        size = ksize;
      } else {
        std::cout << "Axis out of bounds" << std::endl;
        return(0);
      }

      // Make geom (clean me)
      Array<Box> geomProbDomain = amrDataPtrArray[iPlot]->ProbDomain();
      Array<Real> probLo=amrDataPtrArray[iPlot]->ProbLo();
      Array<Real> probHi=amrDataPtrArray[iPlot]->ProbHi();
      RealBox rb(probLo.dataPtr(),probHi.dataPtr());
      Array<Geometry> geom(nLevels);
      for (int iLevel=0; iLevel<nLevels; iLevel++)
	geom[iLevel].define(geomProbDomain[iLevel], &rb, 0);

      //
      // Determine the intesects
      //
      for (int iLevel=0; iLevel<finestLevel; iLevel++) {
	  BoxArray baf = mf[iLevel+1].boxArray();
	  baf.coarsen(amrDataPtrArray[iPlot]->RefRatio()[iLevel]);

	  for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {
	      FArrayBox& myFab = mf[iLevel][mfi];

	      myFab.setVal(1.0,nVars); // Sets the intesect flag to one // To be set to zero where finer level data exists

	      int idx = mfi.index();

	      std::vector< std::pair<int,Box> > isects = baf.intersections(mf[iLevel].boxArray()[idx]);

	      for (int ii = 0; ii < isects.size(); ii++) {
		  myFab.setVal(0.0,isects[ii].second,nVars,1);
	      }
	  }
	  mf[iLevel].FillBoundary(0,mf[iLevel].nComp());
	  mf[iLevel].EnforcePeriodicity(0,nVars,geom[iLevel].periodicity());
      }

      for (MFIter mfi(mf[finestLevel]); mfi.isValid(); ++mfi)
	  mf[finestLevel][mfi].setVal(1.0,nVars); // Sets the intesect flag to one for all fine data

      //
      // Now we have the multi-level data, with an intersect flag
      // Let's integrate...
      //

      int refRatio = 1;

      // Hack next two lines to taste - need to do levTurb below too
      int turbVars(nVars); // let's default to integrating each var coming in
      Real *turb = (Real*)malloc(turbVars*size*sizeof(Real));

      if (turb==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
      for (int i=0; i<turbVars*size; i++) turb[i]=0.;

      // Let's go...
      if (ParallelDescriptor::IOProcessor())
	  //std::cout << "Integrating..." << std::endl;
      for (int iLevel=finestLevel; iLevel>=0; iLevel--) {
	  //if (ParallelDescriptor::IOProcessor())
	      //std::cout << "   Level = " << iLevel << std::endl;
	  Box levProbDomain = amrDataPtrArray[iPlot]->ProbDomain()[iLevel];
	  int levIsize(levProbDomain.length(Amrvis::XDIR));
	  int levJsize(levProbDomain.length(Amrvis::YDIR));
	  int levKsize(levProbDomain.length(Amrvis::ZDIR));
	  int levsize;
	  if (axis == 0) {
	    levsize = levIsize;
	  } else if (axis == 1) {
	    levsize = levJsize;
	  } else if (axis == 2) {
	    levsize = levKsize;
	  } else {
	    std::cout << "Axis out of bounds" << std::endl;
	    return(0);
	  }
	  // grab some ram
	  Real *levTurb = (Real*)malloc(turbVars*levsize*sizeof(Real));

	  if (levTurb==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
	  for (int i=0; i<turbVars*levsize; i++) levTurb[i]=0.;

	  for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {

	      FArrayBox& myFab = mf[iLevel][mfi];

	      int idx = mfi.index();

	      const Real* dataPtr = myFab.dataPtr();
	      const int*  dlo     = myFab.loVect();
	      const int*  dhi     = myFab.hiVect();
	      const int*  lo      = mf[iLevel].boxArray()[idx].loVect();
	      const int*  hi      = mf[iLevel].boxArray()[idx].hiVect();
	      const Real* levDx   = amrDataPtrArray[iPlot]->DxLevel()[iLevel].dataPtr();

	      FORT_INTEGRATE(dataPtr,&nVars,
			     ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
			     levDx,levTurb,&levIsize,&levJsize,&levKsize,&turbVars,&axis,&size);
	  }

	  if (iLevel<finestLevel)  refRatio *= amrDataPtrArray[iPlot]->RefRatio()[iLevel];

	  if (ParallelDescriptor::IOProcessor())
	      //std::cout << "      refRatio = " << refRatio << std::endl;

	  // map levTurb onto turb
	  for (int l=0, k=0; l<levsize; l++)
	      for (int r=0; r<refRatio; r++, k++)
		  for (int v=0; v<turbVars; v++)
		      turb[k*turbVars+v] += levTurb[l*turbVars+v];

	  free(levTurb);
      }

      ParallelDescriptor::ReduceRealSum(turb, size*turbVars, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
	  // Divide by the area
	  Real Lx = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::XDIR];
 	  Real Ly = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::YDIR];
 	  Real Lz = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::ZDIR];
	  Real area;
	  Real vol;

	  if (axis == 0) {
	    area = Ly*Lz;
	    vol = area*Lx;
	  } else if (axis == 1) {
	    area = Lx*Lz;
	    vol = area*Ly;
	  } else if (axis == 2) {
	    area = Lx*Ly;
	    vol = area*Lz;
	  } else {
	    std::cout << "Axis out of bounds" << std::endl;
	    return(0);
	  }

	  for (int i=0; i<turbVars*size; i++) {
	    turb[i] /= area;
	  }

	  Real *axisArr = (Real*)malloc(size*sizeof(Real));
	  for (int i=0; i<size; i++) {
	      axisArr[i] = probLo[axis]+dx[axis]*(0.5+(Real)i);
	  }




	  for (int v=0; v<turbVars; v++) {
	    Real integral = 0.;
	    for (int k=0; k<size; k++) {
	      //std::cout << turb[k*turbVars+v] << std::endl;
	      //integral += turb[k*turbVars+v];
	    }
	    //std::cout << "Integral(" << v << ") = " << integral*vol/(Real)size << std::endl;
	  }
	  //
	  // Write data file
	  //
	  FILE *file;
	  char filename[512];
	  sprintf(filename,"%s/intMean.dat",plotFileNames[iPlot].c_str());
	  file = fopen(filename,"w");
	  fwrite(&(time[iPlot]),sizeof(Real),1,file);
	  fwrite(       &size ,sizeof(int),1,file);
	  fwrite(    &turbVars ,sizeof(int),1,file);
	  fwrite(            axisArr ,sizeof(Real),size,file);
	  fwrite(         turb ,sizeof(Real),size*turbVars,file);
	  fclose(file);

	  free(axisArr);
    //std::cout << "Finish mean calc" << std::endl;
      }



      //Standard deviation calcs
      refRatio = 1;
      Real *turbStd = (Real*)malloc(turbVars*size*sizeof(Real));

      if (turbStd==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
      for (int i=0; i<turbVars*size; i++) turbStd[i]=0.;

       for (int iLevel=finestLevel; iLevel>=0; iLevel--) {
	  //if (ParallelDescriptor::IOProcessor())
	      //std::cout << "   Level = " << iLevel << std::endl;
	  Box levProbDomain = amrDataPtrArray[iPlot]->ProbDomain()[iLevel];
	  int levIsize(levProbDomain.length(Amrvis::XDIR));
	  int levJsize(levProbDomain.length(Amrvis::YDIR));
	  int levKsize(levProbDomain.length(Amrvis::ZDIR));
	  int levsize;
	  if (axis == 0) {
	    levsize = levIsize;
	  } else if (axis == 1) {
	    levsize = levJsize;
	  } else if (axis == 2) {
	    levsize = levKsize;
	  } else {
	    std::cout << "Axis out of bounds" << std::endl;
	    return(0);
	  }
	  // grab some ram
	  Real *levTurb = (Real*)malloc(turbVars*levsize*sizeof(Real));

	  if (levTurb==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
	  for (int i=0; i<turbVars*levsize; i++) levTurb[i]=0.;

	  for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {

	      FArrayBox& myFab = mf[iLevel][mfi];

	      int idx = mfi.index();

	      const Real* dataPtr = myFab.dataPtr();
	      const int*  dlo     = myFab.loVect();
	      const int*  dhi     = myFab.hiVect();
	      const int*  lo      = mf[iLevel].boxArray()[idx].loVect();
	      const int*  hi      = mf[iLevel].boxArray()[idx].hiVect();
	      const Real* levDx   = amrDataPtrArray[iPlot]->DxLevel()[iLevel].dataPtr();

	      FORT_STDDEV(dataPtr,&nVars,
			     ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
			     levDx,levTurb,&levIsize,&levJsize,&levKsize,&turbVars,&axis,&size,turb);
	  }

	  if (iLevel<finestLevel)  refRatio *= amrDataPtrArray[iPlot]->RefRatio()[iLevel];

	  if (ParallelDescriptor::IOProcessor())
	      //std::cout << "      refRatio = " << refRatio << std::endl;

	  // map levTurb onto turb

	  for (int l=0, k=0; l<levsize; l++)
	      for (int r=0; r<refRatio; r++, k++)
		  for (int v=0; v<turbVars; v++)
		      turbStd[k*turbVars+v] += levTurb[l*turbVars+v];

	  free(levTurb);
      }
      ParallelDescriptor::ReduceRealSum(turbStd, size*turbVars, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
	  // Divide by the number of points
	  int numpoints;

	  if (axis == 0) {
	    numpoints = jsize*ksize;
	  } else if (axis == 1) {
	    numpoints = isize*ksize;
	  } else if (axis == 2) {
	    numpoints = isize*jsize;
	  } else {
	    std::cout << "Axis out of bounds" << std::endl;
	    return(0);
	  }

	  for (int i=0; i<turbVars*size; i++) {
	    turbStd[i] /= numpoints;
      turbStd[i] = sqrt(turbStd[i]);
	  }

	  Real *axisArr = (Real*)malloc(size*sizeof(Real));
	  for (int i=0; i<size; i++) {
	      axisArr[i] = probLo[axis]+dx[axis]*(0.5+(Real)i);
	  }




    for (int k=0; k<size; k++) {
	    Real integral = 0.;
      std::cout << axisArr[k] << " ";
	    for (int v=0; v<turbVars; v++) {
	      std::cout << turb[k*turbVars+v] << " ";
	      //integral += turb[k*turbVars+v];
	    }
      std::cout << std::endl;
	    //std::cout << "Integral(" << v << ") = " << integral*vol/(Real)size << std::endl;
	  }
	  //
	  // Write data file
	  //
	  FILE *file;
	  char filename[512];
	  sprintf(filename,"%s/stdDev.bin",plotFileNames[iPlot].c_str());
	  file = fopen(filename,"w");
	  fwrite(&(time[iPlot]),sizeof(Real),1,file);
	  fwrite(       &size ,sizeof(int),1,file);
	  fwrite(    &turbVars ,sizeof(int),1,file);
	  fwrite(            axisArr ,sizeof(Real),size,file);
	  fwrite(         turbStd ,sizeof(Real),size*turbVars,file);
	  fclose(file);

	  free(axisArr);
      }










      free(turb);
      free(turbStd);
      ParallelDescriptor::Barrier();
  }

  BoxLib::Finalize();
  return 0;
}

static
void
PrintUsage (char* progName)
{
  std::cout << "\nUsage:\n"
	    << progName
	    << "\n\tinfile = inputFileName"
	    << "\n\t[-help]"
	    << "\n\n";
  exit(1);
}
