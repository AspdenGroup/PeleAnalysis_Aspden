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

    int nVars(3);
    Array<string> whichVar(nVars);
    whichVar[0] = "x_velocity";
    whichVar[1] = "y_velocity";
    whichVar[2] = "z_velocity";

    Array<int>    destFills(nVars);
    for (int c=0; c<nVars; c++)
	destFills[c] = c;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
    Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
    Array<Real>           time(nPlotFiles);
    Array<Real>           urms(nPlotFiles);

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
	Real dxyz = dx[0]*dx[1]*dx[2];

	int ngrow(0);
	MultiFab mf;
	mf.define(amrDataPtrArray[iPlot]->boxArray(finestLevel), nVars, ngrow, Fab_allocate);

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "Processing " << iPlot << "/" << nPlotFiles << std::endl;
	amrDataPtrArray[iPlot]->FillVar(mf, finestLevel, whichVar, destFills);
	for (int n=0; n<nVars; n++)
	    amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(whichVar[n]));

	Real vol(0), uxb(0), uyb(0), uzb(0), ux2(0), uy2(0), uz2(0);
	for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
	    const FArrayBox &myFab = mf[ntmfi];
	    Array<const Real *> varPtr(nVars);
	    for (int v=0; v<nVars; v++)
		varPtr[v] = myFab.dataPtr(v);
	    const int  *lo  = ntmfi.validbox().loVect();
	    const int  *hi  = ntmfi.validbox().hiVect();

	    int ix = hi[0]-lo[0]+1;
	    int jx = hi[1]-lo[1]+1;
	    int kx = hi[2]-lo[2]+1;
	    for (int k=0; k<kx; k++) {
		Real z=probLo[2] + dx[2]*(0.5+(Real)(k+lo[2]));
		for (int j=0; j<jx; j++) {
		    Real y=probLo[1] + dx[1]*(0.5+(Real)(j+lo[1]));
		    for (int i=0; i<ix; i++) {
			Real x=probLo[0] + dx[0]*(0.5+(Real)(i+lo[0]));
			int cell = (k*jx+j)*ix+i;
			Real ux = varPtr[0][cell];
			Real uy = varPtr[1][cell];
			Real uz = varPtr[2][cell];
			vol += dxyz;
			uxb += ux*dxyz;
			uyb += uy*dxyz;
			uzb += uz*dxyz;
			ux2 += ux*ux*dxyz;
			uy2 += uy*uy*dxyz;
			uz2 += uz*uz*dxyz;
		    }
		}
	    }
	}
	ParallelDescriptor::ReduceRealSum(vol);
	ParallelDescriptor::ReduceRealSum(uxb);
	ParallelDescriptor::ReduceRealSum(uyb);
	ParallelDescriptor::ReduceRealSum(uzb);
	ParallelDescriptor::ReduceRealSum(ux2);
	ParallelDescriptor::ReduceRealSum(uy2);
	ParallelDescriptor::ReduceRealSum(uz2);
	uxb /= vol;	uyb /= vol;	uzb /= vol;
	ux2 /= vol;	uy2 /= vol;	uz2 /= vol;
	urms[iPlot] = sqrt( ( (ux2-uxb*uxb) + (uy2-uyb*uyb) + (uz2-uzb*uzb) ) / 3. );
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << "   ...done." << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
	FILE *file = fopen("RmsVel.dat","w");
	for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
	    fprintf(file,"%e %e\n",time[iPlot],urms[iPlot]);
	fclose(file);
    }

    BoxLib::Finalize();
    return 0;
}

