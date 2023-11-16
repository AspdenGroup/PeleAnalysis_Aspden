//
// This is a version of AmrDerive.cpp that calculates integrals of
// quantities and writes out scalars instead of plotfiles.
//

#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>

using namespace amrex;

extern "C"
{
  void FORT_INTEGRATE(const Real* dat,
			const int* dlo, const int* dhi,
			const int* lo, const int* hi,
			const Real* xlo, const Real* delta, Real* turb,
			int *rsize, int *levRsize,int *ksize, int *levKsize,
		        int *nVars, int *turbVars);
  void pushvtog(const int* lo,  const int* hi,
                const int* dlo, const int* dhi,
                Real* U, const int* Ulo, const int* Uhi,
                const int* nc);
} 


int main (int   argc, char* argv[])
{
  Initialize(argc,argv);
  ParmParse pp;
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
  Vector<std::string> plotFileNames(nPlotFiles);
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    pp.get("infile", plotFileNames[iPlot], iPlot);
    //VSHOWVAL(verbose, plotFileNames[iPlot]);
  }

  // More random initialisation stuff
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  Vector<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
  Vector<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
  Vector<Real>           time(nPlotFiles);

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
  
  //int nVars(pp.countval("vars"));
  int nVars= 5;
  Vector<std::string> varName(nVars);
  varName[0] = "x_velocity";
  varName[1] = "y_velocity";
  varName[2] = "z_velocity";
  varName[3] = "density";
  varName[4] = "temp";
  //varName[5] = "Y(H2)";
  //varName[6] = "Y(O2)";
  //varName[7] = "Y(N2)";
  //varName[5] = "HeatRelease";
  Vector<int> destFillComps(nVars);
  for (int i=0; i<nVars; ++i)
    destFillComps[i] = i;

  // only necessary to set to 1 if derivatives required
  int nGrow = 1;

  //
  // Loop over plots
  //
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

      //
      // Load the data
      //
      Vector<MultiFab> mf(nLevels);
      for (int iLevel=0; iLevel<nLevels; iLevel++) {
	  const BoxArray& ba = amrDataPtrArray[iPlot]->boxArray(iLevel);
	  mf[iLevel].define(ba,DistributionMapping(ba),nVars+1,nGrow); // Do nVars+1 so that we can use the last one as a flag
	  if (ParallelDescriptor::IOProcessor())
	      std::cout << "Reading level: " << iLevel << std::endl;
	  amrDataPtrArray[iPlot]->FillVar(mf[iLevel],iLevel,varName,destFillComps);   // FillVar
	  if (ParallelDescriptor::IOProcessor())
	      std::cout << "...done reading level: " << iLevel << std::endl;
	  for (int iVar=0; iVar<nVars; iVar++)
	      amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varName[iVar]));   // Flush grids
      }
      
      if (nGrow) {
	for (int iLevel=0; iLevel <= finestLevel; iLevel++) {
	const Box& dbox = amrDataPtrArray[iPlot]->ProbDomain()[iLevel];
	  for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {
	      const Box& box = mfi.validbox();
	      FArrayBox& fab = mf[iLevel][mfi];
	      int ncomp = mf[iLevel].nComp();
	      pushvtog(BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_BOX(dbox),
               BL_TO_FORTRAN_ANYD(fab),
               &ncomp);
	    }	 
	  mf[iLevel].FillBoundary();
	}
      }
		 
      
      const Real *dx = amrDataPtrArray[iPlot]->DxLevel()[finestLevel].dataPtr();
      Box probDomain = amrDataPtrArray[iPlot]->ProbDomain()[finestLevel];

      int isize(probDomain.length(Amrvis::XDIR));
      int jsize(probDomain.length(Amrvis::YDIR));
#if (BL_SPACEDIM==3)
      int ksize(probDomain.length(Amrvis::ZDIR));
#endif
      //std::cout << "ksize = " << ksize << std::endl;
      int rsize = (int)(std::min(isize,jsize)/2.0);
      //std::cout << "rsize = " << rsize << std::endl;
      
      // Make geom (clean me)
      Vector<Box> geomProbDomain = amrDataPtrArray[iPlot]->ProbDomain();
      Vector<Real> probLo=amrDataPtrArray[iPlot]->ProbLo();
      Vector<Real> probHi=amrDataPtrArray[iPlot]->ProbHi();
      RealBox rb(probLo.dataPtr(),probHi.dataPtr());
      Vector<Geometry> geom(nLevels);
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
      turbVars = 6; //Actually this many now
      Real *turb = (Real*)malloc(turbVars*ksize*rsize*sizeof(Real));

      if (turb==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
      for (int i=0; i<turbVars*ksize*rsize; i++) turb[i]=0.;

      // Let's go...
      if (ParallelDescriptor::IOProcessor())
	  std::cout << "Integrating..." << std::endl;
      for (int iLevel=finestLevel; iLevel>=0; iLevel--) {
	  if (ParallelDescriptor::IOProcessor())
	     std::cout << "   Level = " << iLevel << std::endl;
	  Box levProbDomain = amrDataPtrArray[iPlot]->ProbDomain()[iLevel];
	  int levIsize(levProbDomain.length(Amrvis::XDIR));
	  int levJsize(levProbDomain.length(Amrvis::YDIR));
#if (BL_SPACEDIM==3)
	  int levKsize(levProbDomain.length(Amrvis::ZDIR));
#endif
	  int levRsize = std::min({levIsize,levJsize})/2;
	  //int levRsize = (int)(std::sqrt(levIsize*levIsize + levJsize*levJsize)/2.0);
	  // grab some ram
	  /* APPARENTLY FORTRAN DOES THIS NOW
	  Real *levTurb = (Real*)malloc(turbVars*levrsize*sizeof(Real));

	  if (levTurb==NULL) {std::cout << "Error: Couldn't allocate turb stats memory" << std::endl; return(0);}
	  for (int i=0; i<turbVars*levsize; i++) levTurb[i]=0.;
	  */
	  for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {

	      FArrayBox& myFab = mf[iLevel][mfi];

	      int idx = mfi.index();
	      const Real* dataPtr = myFab.dataPtr();
	      const int*  dlo     = myFab.smallEnd().getVect();
	      const int*  dhi     = myFab.bigEnd().getVect();
	      const Box& box = mfi.validbox(); 
	      const int*  lo      = box.smallEnd().getVect();
	      const int*  hi      = box.bigEnd().getVect();
	      const Real* xlo     = amrDataPtrArray[iPlot]->GridLocLo()[iLevel][idx].dataPtr();
	      const Real* levDx   = amrDataPtrArray[iPlot]->DxLevel()[iLevel].dataPtr();

	      FORT_INTEGRATE(dataPtr,
			     dlo,dhi,lo,hi,xlo,
			     levDx,turb,&rsize,&levRsize,&ksize,&levKsize,&nVars,&turbVars);
	  }
	  //if (iLevel<finestLevel)  refRatio *= amrDataPtrArray[iPlot]->RefRatio()[iLevel];
	  /* DONE IN FORTRAN
	  if (ParallelDescriptor::IOProcessor())
	    std::cout << "      refRatio = " << refRatio << std::endl;

	  // map levTurb onto turb
	  
	  for (int l=0, k=0; l<levsize; l++)
	    for (int r=0; r<refRatio; r++, k++)
	      for (int v=0; v<turbVars; v++)
		turb[k*turbVars+v] += levTurb[l*turbVars+v];
	  
	  free(levTurb);
	  */
      }


      ParallelDescriptor::ReduceRealSum(turb, ksize*rsize*turbVars, ParallelDescriptor::IOProcessorNumber());
      
      if (ParallelDescriptor::IOProcessor()) {
	  // Divide by the area
	  Real Lx = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::XDIR];
 	  Real Ly = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::YDIR];
#if (BL_SPACEDIM==3)
	  Real Lz = amrDataPtrArray[iPlot]->ProbSize()[Amrvis::ZDIR];
#endif
	  
	  //Real area = Lz*std::min({Lx,Ly}); //Need to discuss this
	  

	  for (int i=0; i<turbVars*ksize*rsize; i++) {
	    //turb[i] /= area;//*turb[i-i%turbVars];
	    if (i%turbVars != 0) {
	      turb[i] /= turb[i-i%turbVars];
	    }
	   }
	  

	  Real *r = (Real*)malloc(rsize*sizeof(Real));
	  for (int i=0; i<rsize; i++) {
	    r[i] = std::sqrt(2)*dx[0]*(0.5+(Real)i);
	  }
	  
	  Real *z = (Real*)malloc(ksize*sizeof(Real));
	  for (int k=0; k<ksize; k++) {
	      z[k] = dx[2]*(0.5+(Real)k);
	  }
	  //
	  // Write data file
	  //
	  FILE *file;
	  char filename[512];
	  for (int k = 0; k < ksize; k++) {
	    //fprintf(file,"%e ",z[k]);
	    sprintf(filename,"%s/radMeanz%d.dat",plotFileNames[iPlot].c_str(),k);
	    file = fopen(filename,"w");
	    for (int j = 0; j < rsize; j++) {
	      fprintf(file,"%e ",r[j]);
	      for (int v=0; v<turbVars;v++) {
		fprintf(file,"%e ",turb[(k*rsize+j)*turbVars+v]);
	      }
	      fprintf(file,"\n");
	    }
	    fclose(file);
          } 
	  free(r);
	  free(z);
      }
      free(turb);

      ParallelDescriptor::Barrier();
      }
	

  Finalize();
  return 0;
  }


  static void
    PrintUsage (char* progName)
  {
  std::cout << "\nUsage:\n"
	    << progName
	    << "\n\tinfile = inputFileName"
	    << "\n\t[-help]"
	    << "\n\n";
  exit(1);
  }

