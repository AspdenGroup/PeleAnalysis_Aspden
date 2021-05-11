
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

#include "AmrDeriveGenericIntegration_F.H"

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
  int nProcs(ParallelDescriptor::NProcs());    VSHOWVAL(verbose, nProcs);
  int myProc = ParallelDescriptor::MyProc();

  // Count number of input plot files
  int nPlotFiles(pp.countval("infile"));
  if(nPlotFiles <= 0) {
    std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
    std::cerr << "Exiting." << std::endl;
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  VSHOWVAL(verbose, nPlotFiles);

  // Make an array of srings containing paths of input plot files
  Array<std::string> plotFileNames(nPlotFiles);
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    pp.get("infile", plotFileNames[iPlot], iPlot);
    VSHOWVAL(verbose, plotFileNames[iPlot]);
  }

  // More random initialisation stuff
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
  Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
  Array<Real>           time(nPlotFiles);
  Array<Real>           int_volume(nPlotFiles);
    
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Loading " << plotFileNames[iPlot] << std::endl;
      
    dataServicesPtrArray[iPlot] = new DataServices(plotFileNames[iPlot], fileType);               // Populate DataServices array
      
    if( ! dataServicesPtrArray[iPlot]->AmrDataOk())                                               // Check AmrData ok
      DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not
      
    amrDataPtrArray[iPlot] = &(dataServicesPtrArray[iPlot]->AmrDataRef());                        // Populate DataPtrArray

    time[iPlot] = amrDataPtrArray[iPlot]->Time();

    if (verbose) std::cout << "Time = " << time[iPlot] << std::endl;
  }
    
  int finestLevel = amrDataPtrArray[0]->FinestLevel();
  pp.query("finestLevel",finestLevel);
  int nLevels = finestLevel + 1;

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
	      std::cout << "Reading level: " << iLevel << std::endl;
	  amrDataPtrArray[iPlot]->FillVar(mf[iLevel],iLevel,varName,destFillComps);   // FillVar
	  if (ParallelDescriptor::IOProcessor())
	      std::cout << "...done reading level: " << iLevel << std::endl;
	  for (int iVar=0; iVar<nVars; iVar++)
	      amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varName[iVar]));   // Flush grids
      }

      const Real *dx = amrDataPtrArray[iPlot]->DxLevel()[finestLevel].dataPtr();
      Box probDomain = amrDataPtrArray[iPlot]->ProbDomain()[finestLevel];

      int isize(probDomain.length(Amrvis::XDIR));
      int jsize(probDomain.length(Amrvis::YDIR));
#if (BL_SPACEDIM==2)
      int size = jsize;
#else
      int ksize(probDomain.length(Amrvis::ZDIR));
      int size = ksize;
#endif
      
      
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
	  std::cout << "Integrating..." << std::endl;
      for (int iLevel=finestLevel; iLevel>=0; iLevel--) {
	  if (ParallelDescriptor::IOProcessor())
	      std::cout << "   Level = " << iLevel << std::endl;
	  Box levProbDomain = amrDataPtrArray[iPlot]->ProbDomain()[iLevel];
	  int levIsize(levProbDomain.length(Amrvis::XDIR));
	  int levJsize(levProbDomain.length(Amrvis::YDIR));
#if (BL_SPACEDIM==2)
	  int levsize = levJsize;
#else
	  int levKsize(levProbDomain.length(Amrvis::ZDIR));
	  int levsize = levKsize;
#endif
	  
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

#if (BL_SPACEDIM==2)
	      FORT_INTEGRATE(dataPtr,&nVars,
			     ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
			     levDx,levTurb,&levIsize,&levJsize,&turbVars);
#else
	      FORT_INTEGRATE(dataPtr,&nVars,
			     ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
			     levDx,levTurb,&levIsize,&levJsize,&levKsize,&turbVars);
#endif
	  }
	  
	  if (iLevel<finestLevel)  refRatio *= amrDataPtrArray[iPlot]->RefRatio()[iLevel];

	  if (ParallelDescriptor::IOProcessor())
	      std::cout << "      refRatio = " << refRatio << std::endl;
	  
	  // map levTurb onto turb
#if (BL_SPACEDIM==2)
	  for (int l=0, k=0; l<levJsize; l++)
	      for (int r=0; r<refRatio; r++, k++)
		  for (int v=0; v<turbVars; v++)
		      turb[k*turbVars+v] += levTurb[l*turbVars+v];
#else
	  for (int l=0, k=0; l<levKsize; l++)
	      for (int r=0; r<refRatio; r++, k++)
		  for (int v=0; v<turbVars; v++)
		      turb[k*turbVars+v] += levTurb[l*turbVars+v];
#endif	  
	  free(levTurb);
      }

      ParallelDescriptor::ReduceRealSum(turb, size*turbVars, ParallelDescriptor::IOProcessorNumber());
      
      if (ParallelDescriptor::IOProcessor()) {
	  // Divide by the area
	  Real Lx = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::XDIR];
 	  Real Ly = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::YDIR];
#if (BL_SPACEDIM==2)
	  Real area = Lx;
	  Real vol = area * Ly;
#else
 	  Real Lz = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::ZDIR];
	  Real area = Lx * Ly;
	  Real vol = area*Lz;
#endif
	  for (int i=0; i<turbVars*size; i++) {
	    turb[i] /= area;
	  }

	  for (int v=0; v<turbVars; v++) {
	    Real integral = 0.;
#if (BL_SPACEDIM==2)
	    for (int j=0; j<jsize; j++)
	      integral += turb[j*turbVars+v];
	    std::cout << "Integral(" << v << ") = " << integral*vol/(Real)jsize << std::endl;
	    int_volume[iPlot] = integral*vol/(Real)jsize;
#else
	    for (int k=0; k<ksize; k++)
	      integral += turb[k*turbVars+v];
	    std::cout << "Integral(" << v << ") = " << integral*vol/(Real)ksize << std::endl;
	    int_volume[iPlot] = integral*vol/(Real)ksize;
#endif
	  }

	  

	  Real *z = (Real*)malloc(size*sizeof(Real));
	  for (int i=0; i<size; i++)
	      z[i] = dx[2]*(0.5+(Real)i);

 	  //
	  // Write data file
	  //
	  if (ParallelDescriptor::IOProcessor())
	    {
	      FILE *file = fopen("int.dat","w");
	      for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
		fprintf(file,"%e %e\n",time[iPlot],int_volume[iPlot]);
	      fclose(file);
	    }
	  
	  free(z);
      }
      
      free(turb);

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


