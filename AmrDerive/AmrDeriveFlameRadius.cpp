#include <winstd.H>

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <fstream>

#ifndef WIN32
#include <unistd.h>
#endif

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Utility.H"
#include "ArrayLim.H"

static
void
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<filename>\n";
    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);
    
    int nPlotFiles(pp.countval("infile"));
    Array<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot)
	pp.get("infile", plotFileNames[iPlot], iPlot);

    int nVars(1);
    Array<string> whichVar(nVars);
    whichVar[0] = "Y(H2O)";

    Array<int>    destFills(nVars);
    for (int c=0; c<nVars; c++)
	destFills[c] = c;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
    Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
    Array<Real>           time(nPlotFiles);
    Array<Real>           eff_rad(nPlotFiles);

    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Loading " << plotFileNames[iPlot] << std::endl;

	dataServicesPtrArray[iPlot] = new DataServices(plotFileNames[iPlot], fileType);               // Populate DataServices array

	if( ! dataServicesPtrArray[iPlot]->AmrDataOk())                                               // Check AmrData ok
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not

	amrDataPtrArray[iPlot] = &(dataServicesPtrArray[iPlot]->AmrDataRef());                        // Populate DataPtrArray

	time[iPlot] = amrDataPtrArray[iPlot]->Time();

    }

    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

	int finestLevel = amrDataPtrArray[iPlot]->FinestLevel();
	int inFinestLevel(-1);    pp.query("finestLevel",inFinestLevel);
	if (inFinestLevel>-1 && inFinestLevel<finestLevel) {
	    finestLevel = inFinestLevel;
            if (ParallelDescriptor::IOProcessor())
	        std::cout << "Finest level: " << finestLevel << std::endl;
	}

	Array<Real> probLo=amrDataPtrArray[iPlot]->ProbLo();
	Array<Real> probHi=amrDataPtrArray[iPlot]->ProbHi();
	const Real *dx = amrDataPtrArray[iPlot]->DxLevel()[finestLevel].dataPtr();
	Real dxyz;
	int ngrow(0);
	MultiFab mf;
	mf.define(amrDataPtrArray[iPlot]->boxArray(finestLevel), nVars, ngrow, Fab_allocate);

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "Processing " << iPlot << "/" << nPlotFiles << std::endl;
	amrDataPtrArray[iPlot]->FillVar(mf, finestLevel, whichVar, destFills);
	for (int n=0; n<nVars; n++)
	    amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(whichVar[n]));

	Real vol(0);
	for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
	    const FArrayBox &myFab = mf[ntmfi];
	    Array<const Real *> varPtr(nVars);
	    for (int v=0; v<nVars; v++)
		varPtr[v] = myFab.dataPtr(v);
	    const int  *lo  = ntmfi.validbox().loVect();
	    const int  *hi  = ntmfi.validbox().hiVect();
	    int ix = hi[0]-lo[0]+1;
	    int jx = hi[1]-lo[1]+1;
#if (BL_SPACEDIM==3)	      
	    int kx = hi[2]-lo[2]+1;
#endif

#if (BL_SPACEDIM==3)
	    for (int k=0; k<kx; k++) {
		Real z=probLo[2] + dx[2]*(0.5+(Real)(k+lo[2]));
		for (int j=0; j<jx; j++) {
		    Real y=probLo[1] + dx[1]*(0.5+(Real)(j+lo[1]));
		    for (int i=0; i<ix; i++) {
			Real x=probLo[0] + dx[0]*(0.5+(Real)(i+lo[0]));
			int cell = (k*jx+j)*ix+i;
			Real H2O_val = varPtr[0][cell];
			if ( H2O_val > 0.01)
			  {
			    dxyz = dx[0]*dx[1]*dx[2];
			    std::cout << "Cell Volume =" << dxyz << std::endl;
			  }
			else
			  {
			    dxyz = 0;
			    std::cout << "Cell Volume =" << dxyz << std::endl;
			  }

			vol += dxyz;
		    }
		}
	    }
#endif
#if (BL_SPACEDIM==2)

	    for (int j=0; j<jx; j++) {
	      Real y=probLo[1] + dx[1]*(0.5+(Real)(j+lo[1]));
	      for (int i=0; i<ix; i++) {
		Real x=probLo[0] + dx[0]*(0.5+(Real)(i+lo[0]));
		int cell = j*ix+i;
		Real H2O_val = varPtr[0][cell];
		
		if (H2O_val > 0.01)
		  {
		    
		    dxyz = dx[0] * dx[1];
		    //std::cout << "Cell Volume = " << dxyz << std::endl;
		    
		  }
		else
		  {
		    dxyz = 0;
		  }
		vol += dxyz;
	      }
	    }
#endif
	    
	    
	}
	std::cout << "Flame Volume = " << vol << std::endl;
	ParallelDescriptor::ReduceRealSum(vol);
#if (BL_SPACEDIM==2)
	eff_rad[iPlot] = pow((vol / M_PI),0.5);
#endif
#if (BL_SPACEDIM==3)
	eff_rad[iPlot] = pow(( (3/4) * (vol / M_PI) ),(1/3));
#endif
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << "   ...done." << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
	FILE *file = fopen("FlameRadius.dat","w");
	for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
	    fprintf(file,"%e %e\n",time[iPlot],eff_rad[iPlot]);
	fclose(file);
    }

    BoxLib::Finalize();
    return 0;
}
