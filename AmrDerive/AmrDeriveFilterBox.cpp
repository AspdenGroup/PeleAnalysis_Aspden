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
#include "Geometry.H"

#include "AmrDeriveFilter_F.H"

// This one is for the ghost cells (filcc)
#include "xtra_F.H"

static
void
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << "\n"
              << "   fSize= val           : Filter size (total filter width will be 2*val+1)\n"
              << "   finestLevel=n        : Finest level\n"
              << "   fuelName= species    : Fuel name (CH4)\n"
              << "   outSuffix=str        : Suffix to add to the plfile name as an alt dir for results (default=\"_fSize%02i\")\n"
              << "   infile= plt1 plt2    : List of plot files" << std::endl;
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

    // Get plot file
    int nPlotFiles(pp.countval("infile"));
    if(nPlotFiles <= 0) {
        std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
        std::cerr << "Exiting." << std::endl;
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // Get filter size
    int fSize(3);
    pp.query("fSize",fSize);

    // Make an array of srings containing paths of input plot files
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Processing " << nPlotFiles << " plotfiles..." << std::endl;
    Array<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
        pp.get("infile", plotFileNames[iPlot], iPlot);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "   " << plotFileNames[iPlot] << std::endl;
    }

    // Get finest level argument
    int inFinestLevel(-1);
    pp.query("finestLevel",inFinestLevel);

    // Note hacked to finestLevel 0
    int iLevel0=0;
    int nGrow0=0;

    // Variables to load
    int nVars(pp.countval("vars"));
    Array<string> whichVar(nVars);
    for(int v = 0; v < nVars; v++) {
	pp.get("vars", whichVar[v], v);
    }
    Array<int> destFills(nVars);
    for (int c=0; c<nVars; c++)
	destFills[c] = c;

    // Do Favre averaging
    int doFavre(1);
    pp.query("doFavre",doFavre);
    Array<int> favre(nVars);
    for(int v = 0; v < nVars; v++) {
      favre[v]=1;
      if ( (doFavre==0)
	   || (whichVar[v]=="density")
	   || (whichVar[v]=="temp")
	   || (whichVar[v].find("ConsumptionRate")!=std::string::npos)
	   || (whichVar[v].find("HoReact")!=std::string::npos)
	   || (whichVar[v].find("TotHofR")!=std::string::npos)
	   || (whichVar[v]=="HeatRelease")
	   || (whichVar[v]=="phi") )
	favre [v]=0;
    }

    // Out suffix
    char fSizeStr[16];
    if (doFavre)
      sprintf(fSizeStr,"_fSize%02iFavre",fSize);
    else
      sprintf(fSizeStr,"_boxSize%02i",fSize);
    std::string outSuffix(fSizeStr);
    pp.query("outSuffix",outSuffix);

    int denVar=-1;
    if (doFavre) {
      for(int v = 0; v < nVars; v++)
	if (whichVar[v]=="density") denVar=v;
      if (denVar<0)
	BoxLib::Error("Couldn't find density for Favre averaging");
    }

    // Check to screen
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Variable list: (favre flag)" << std::endl;
      for(int v = 0; v < nVars; v++) {
	std::cout << "   " << whichVar[v] << "(" << favre[v] << ")" <<std::endl;
      }
    }

    // Loop over files
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

        // Open file and get an amrData pointer
        std::string infile = plotFileNames[iPlot];
        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nOpening " << infile << "..." << std::endl;
        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(infile, fileType);
        if (!dataServices.AmrDataOk())
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        AmrData& amrData = dataServices.AmrDataRef();
        if (ParallelDescriptor::IOProcessor())
            std::cout << "   ...done." << std::endl;

        // Check the names of the variables are present in the plotfile
        Array<int> destFills(nVars);
        for (int v=0; v<nVars; v++) {
            destFills[v] = v;
            if (amrData.StateNumber(whichVar[v])<0) {
                std::string message="Bad variable name ("+whichVar[v]+")";
                BoxLib::Error(message.c_str());
            }
        }

        // Get finest level
        int finestLevel = amrData.FinestLevel();
        if (inFinestLevel>-1 && inFinestLevel<finestLevel) {
            finestLevel = inFinestLevel;
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Finest level: " << finestLevel << std::endl;
        }

	// Hack
	if (ParallelDescriptor::IOProcessor())
	  std::cout << "WARNING: Overriding finest level to zero..." << std::endl;
	finestLevel=0;

        int nLevels = finestLevel+1;

        // Get domain size
        Array<Real> probLo=amrData.ProbLo();
        Array<Real> probHi=amrData.ProbHi();

        // Get probDomain
        Array<Box> probDomain = amrData.ProbDomain();

	// Make geometry (for fill periodic boundary conditions)
	RealBox rb(probLo.dataPtr(),probHi.dataPtr());
	Array<Geometry> geom(nLevels);
	for (int iLevel=0; iLevel<nLevels; iLevel++)
	  geom[iLevel].define(probDomain[iLevel], &rb, 0);

	// Let's check we're getting periodicity right
	if (ParallelDescriptor::IOProcessor()) {
	  std::cout << "isPeridoc = "
		    << geom[iLevel0].isPeriodic(0) << " "
		    << geom[iLevel0].isPeriodic(1) << " "
		    << geom[iLevel0].isPeriodic(2) << " "
		    << std::endl;
	  if ((geom[iLevel0].isPeriodic(0)+geom[iLevel0].isPeriodic(1)+geom[iLevel0].isPeriodic(2))==0)
	    std::cout << "WARNING: no periodic directions - either that's correct or geom not provided" << std::endl;
	}

	//
        // Load the data
	//
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Loading..." << std::endl;

	Real timer_start = ParallelDescriptor::second();

        PArray<MultiFab> mf(nLevels,PArrayManage);

        for (int iLevel=0; iLevel<nLevels; iLevel++) {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "   Level " << iLevel << "..." << std::endl;
	    // Set ngrow to be the filter size
	    int ngrow(fSize);

            // Make space for each component that will be filtered
            mf.set(iLevel, new MultiFab(amrData.boxArray(iLevel), nVars, ngrow));

	    // Load the data
            amrData.FillVar(mf[iLevel], iLevel, whichVar, destFills);

	    // Clean up
	    for (int iVar=0; iVar<nVars; iVar++)
		amrData.FlushGrids(amrData.StateNumber(whichVar[iVar]));
        }

	// Fill boundary conditions with something sensible, then fill internally, and then periodic
        for (int iLevel=0; iLevel<nLevels; iLevel++) {
	    for (MFIter mfi(mf[iLevel]); mfi.isValid(); ++mfi) {
		const Box& box = mfi.validbox();
		FORT_PUSHVTOGFO(box.loVect(), box.hiVect(),
				mf[iLevel][mfi].dataPtr(),
				ARLIM(mf[iLevel][mfi].loVect()),
				ARLIM(mf[iLevel][mfi].hiVect()),&nVars); // only need to do this for nVars
	    }
	    mf[iLevel].FillBoundary(0,nVars);
	    mf[iLevel].EnforcePeriodicity(0,nVars,geom[iLevel].periodicity());
	}

	Real timer_stop = ParallelDescriptor::second();

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "      ...done (" << timer_stop - timer_start << " seconds)." << std::endl;

	//
	//  COARSEN
	//



	// REMOVED - steal from parent directory if you want it back



	//
	// FILTER
	//

        if (ParallelDescriptor::IOProcessor())
	  std::cout<< "Filtering..." << std::endl;

	timer_start = ParallelDescriptor::second();

	// Make the multifab to hold the output data
        MultiFab mfOut;
	mfOut.define(amrData.boxArray(iLevel0), nVars, nGrow0, Fab_allocate);
	mfOut.setVal(0);

	for(MFIter mfi(mf[iLevel0]); mfi.isValid(); ++mfi) {

	  int idx = mfi.index();

	  FArrayBox &myFab   = mf[iLevel0][mfi];
	  Real      *myData  = myFab.dataPtr();
	  FArrayBox &outFab  = mfOut[mfi];
	  Real      *outData = outFab.dataPtr();

	  const int*  dlo    = myFab.loVect();
	  const int*  dhi    = myFab.hiVect();
	  const int*  olo    = outFab.loVect();
	  const int*  ohi    = outFab.hiVect();
	  const Real* dx     = amrData.DxLevel()[iLevel0].dataPtr();


	  FORT_FILTER(myData,ARLIM(dlo),ARLIM(dhi),
		      outData,ARLIM(olo),ARLIM(ohi),
		      probLo.dataPtr(), probHi.dataPtr(), dx,
		      &nVars, &fSize, favre.dataPtr(), &denVar );
	}

	timer_stop = ParallelDescriptor::second();

        if (ParallelDescriptor::IOProcessor())
	  std::cout << "      ...done (" << timer_stop - timer_start << " seconds)." << std::endl;

	//
	// Now do all the work for a plot file
	//
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Outputting..." << std::endl;

	timer_start = ParallelDescriptor::second();

	string pltfile;
	pltfile = infile + outSuffix;

	if (ParallelDescriptor::IOProcessor())
	  if (!BoxLib::UtilCreateDirectory(pltfile, 0755))
	    BoxLib::CreateDirectoryFailed(pltfile);
	ParallelDescriptor::Barrier();

	std::string HeaderFileName = pltfile + "/Header";

	static const std::string the_plot_file_type("NavierStokes-V1.1");

	std::ofstream os;

	if (ParallelDescriptor::IOProcessor()) {
	  os.open(HeaderFileName.c_str());

	  int old_prec = os.precision(15);

	  // The plot file type
	  os << the_plot_file_type << '\n';
	  // The number of variables
	  os << nVars << '\n';
	  // The variable names
	  for(int v=0; v<nVars; v++)
	    os << whichVar[v] << std::endl;
	  // The number of space dimensions
	  os << BL_SPACEDIM << '\n';
	  // Time
	  os << amrData.Time() << '\n';
	  // Finest level
	  os << 0 << '\n';
	  // Domain
	  for (int i=0; i<BL_SPACEDIM; i++)
	    os << amrData.ProbLo()[i] << ' ';
	  os << '\n';
	  for (int i=0; i<BL_SPACEDIM; i++)
	    os << amrData.ProbHi()[i] << ' ';
	  os << '\n';
	  // Refinement ratios
	  os << '\n';
	  // Cell sizes
	  os << amrData.ProbDomain()[0] << ' ';
	  os << '\n';
	  // Time steps
	  os << amrData.LevelSteps()[0] << ' ';
	  os << '\n';
	  // dx
	  for (int j=0; j<BL_SPACEDIM; j++)
	    os << amrData.DxLevel()[0][j] << ' ';
	  os << '\n';
	  // CoordSys & bndry
	  os << "0\n0\n";
	}
	//
	// Build the directory to hold the MultiFab
	// The name is relative to the directory containing the Header file.
	//
	static const std::string BaseName = "/Cell";
	std::string Level = "Level_0";
	std::string FullPath = pltfile;
	if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
	  FullPath += '/';
	FullPath += Level;
	if (ParallelDescriptor::IOProcessor())
	  if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
	    BoxLib::CreateDirectoryFailed(FullPath);
	ParallelDescriptor::Barrier();

	if (ParallelDescriptor::IOProcessor()) {
	  std::cout << "      Grids = " << amrData.boxArray(0).size() << std::endl;

	  os << 0 << ' ' << amrData.boxArray(0).size() << ' ' << amrData.Time() << '\n';
	  os << amrData.LevelSteps()[0] << '\n';

	  for (int i = 0; i < amrData.boxArray(0).size(); i++) {
	    for (int n = 0; n < BL_SPACEDIM; n++)
	      os << amrData.GridLocLo()[0][i][n] << ' ' << amrData.GridLocHi()[0][i][n] << '\n';
	  }
	  std::string PathNameInHeader = Level;
	  PathNameInHeader += BaseName;
	  os << PathNameInHeader << '\n';
	}
	//
	// Use the Full pathname when naming the MultiFab.
	//
	std::string TheFullPath = FullPath;
	TheFullPath += BaseName;
	VisMF::Write(mfOut,TheFullPath);

	timer_stop = ParallelDescriptor::second();

        if (ParallelDescriptor::IOProcessor())
	  std::cout << "      ...done (" << timer_stop - timer_start << " seconds)." << std::endl;

	// Do we need to do some clean up here?

    } // iPlot

    BoxLib::Finalize();
    return 0;
}
