#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> outfile=<> [options] \n\tOptions:\n";
    exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
    vector<std::string> tokens = Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    int numFiles = pp.countval("infiles");
    
    Vector<std::string> infiles; infiles.resize(numFiles);
    pp.getarr("infiles",infiles);
    
    std::string outfile = "avgPlot";
    pp.query("outfile",outfile);
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Print() << "Setting up avgPlotfile..." << std::endl;
    // Set up for reading pltfile
    DataServices dataServices0(infiles[0], fileType);

    if (!dataServices0.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData0 = dataServices0.AmrDataRef();

    // read/write all comps
    int nComp = amrData0.NComp();
    Vector<int> comps(nComp);
    for (int i=0; i<nComp; ++i)
      comps[i] = i;

    int finestLevel = amrData0.FinestLevel();
    //finest level to go down to
    pp.query("finestLevel",finestLevel);
 
    Print() << "... averaging to resolution of level " << finestLevel << std::endl;
    // find probDomain at level we want
    Box probDomain0 = amrData0.ProbDomain()[finestLevel];
    int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    BoxArray tmpnewba(probDomain0);
    tmpnewba.maxSize(max_grid_size);

    const BoxArray& newba = tmpnewba;

    Print() << "... BoxArray set" << std::endl;
    Print() << "... number of boxes = " << newba.size() << std::endl;
    //const DistributionMapping& dm(newba);
    int levelSteps; //dummy variable
    
    RealBox rb(&(amrData0.ProbLo()[0]),
               &(amrData0.ProbHi()[0]));
    Vector<int> is_per(BL_SPACEDIM,0);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    int coord = 0;
    Geometry geoms(probDomain0,&rb,coord,&(is_per[0]));
    Print() << "... creating avgMF" << std::endl;
    std::unique_ptr<MultiFab> avgMF;
    avgMF.reset(new MultiFab(newba,DistributionMapping(newba),nComp,0));
    avgMF->setVal(0);
    Print() << "... avgMF created" << std::endl;
    Real time0 = amrData0.Time();
    Vector<std::string> names(nComp);
    for (int i=0; i<nComp; ++i) {
      names[i] = amrData0.PlotVarNames()[comps[i]];
      amrData0.FlushGrids(comps[i]); //do i need to flush if not filled?
    }
    
    Real timeFinal = 0;
    Print() << "avgPlotfile setup completed." << std::endl;
    Print() << "Starting averaging... " << std::endl;
    for (int iFile = 0; iFile<numFiles; iFile++) {
      // Set up for reading pltfile
      DataServices dataServices(infiles[iFile], fileType);
      
      if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
      
      AmrData& amrData = dataServices.AmrDataRef();

      // make space
      std::unique_ptr<MultiFab> fileData;
      fileData.reset(new MultiFab(newba,DistributionMapping(newba),nComp,0));
	//MultiFab fileData(newba,DistributionMapping(newba),nComp,0);
      Print() << "... loading " << infiles[iFile] << std::endl;
      // load data
      amrData.FillVar(*fileData,finestLevel,amrData.PlotVarNames(),comps);
      for (int i = 0; i < nComp; i++) {
	  amrData.FlushGrids(comps[i]);
      }
      Print() << "... " << infiles[iFile] << " loaded" << std::endl;
      Print() << "... adding to averaged MF" << std::endl;
      MultiFab::Saxpy(*avgMF,1.0/numFiles,*fileData,0,0,nComp,0);
      if(iFile == numFiles-1) {
	timeFinal = amrData.Time();
      }
      Print() << "... deleting fileData" << std::endl;
      fileData.reset();	        
    }
    Print() << "Completed averaging." << std::endl;
    Print() << "Writing " << outfile << "..." << std::endl;
    WriteSingleLevelPlotfile(outfile,*avgMF,names,geoms,timeFinal-time0,levelSteps);
    int doStdDev = 0;
    pp.query("doStdDev",doStdDev);
    if(doStdDev) {
      Print() << "Doing standard deviation" << std::endl;
      std::string stdDevOutfile= "stdDevPlotfile";
      pp.query("stdDevOutfile",stdDevOutfile);
      Print() << "stdDevPlotfile setup..." << std::endl;
      Print() << "... creating stdDevMF" << std::endl;
      std::unique_ptr<MultiFab> stdDevMF;
      stdDevMF.reset(new MultiFab(newba,DistributionMapping(newba),nComp,0));
      stdDevMF->setVal(0);
      Print() << "... stdDevMF created" << std::endl;
      Vector<std::string> names_stddev(nComp);
      for (int i=0; i<nComp; ++i)
	names_stddev[i] = amrData0.PlotVarNames()[comps[i]]+"_stdDev";
      
      Print() << "stdDevPlotfile setup completed." << std::endl;
      Print() << "Starting deviations... " << std::endl;
      for (int iFile = 0; iFile<numFiles; iFile++) {
	// Set up for reading pltfile
	DataServices dataServices(infiles[iFile], fileType);
	
	if (!dataServices.AmrDataOk())
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);
	
	AmrData& amrData = dataServices.AmrDataRef();
	
       // make space
	std::unique_ptr<MultiFab> fileData;
	fileData.reset(new MultiFab(newba,DistributionMapping(newba),nComp,0));
	Print() << "... loading " << infiles[iFile] << std::endl;
	// load data
	amrData.FillVar(*fileData,finestLevel,amrData.PlotVarNames(),comps);
	for (int i = 0; i < nComp; i++) {
	  amrData.FlushGrids(comps[i]);
	}

	Print() << "... " << infiles[iFile] << " loaded" << std::endl;
	Print() << "... subtracting averaged MF" << std::endl;       
	MultiFab::Subtract(*fileData,*avgMF,0,0,nComp,0);
	Print() << "... multiply with self" << std::endl;
	MultiFab::Multiply(*fileData,*fileData,0,0,nComp,0);
	Print() << "... add to stdDevMF" << std::endl;
	MultiFab::Saxpy(*stdDevMF,1.0/numFiles,*fileData,0,0,nComp,0);
	Print() << "... deleting fileData" << std::endl;
	fileData.reset();
      }

      avgMF.reset();
      Print() << "Variation calculated, square rooting..." << std::endl;
      //is this the best way to do this?
      for (MFIter mfi(*stdDevMF,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	Array4<Real> array = stdDevMF->array(mfi);
	AMREX_PARALLEL_FOR_4D (bx, nComp, i, j, k, n,
			       {
				 array(i,j,k,n) = std::sqrt(array(i,j,k,n));
			       });
      }
      Print() << "Finished calculating stardard deviation." << std::endl;
      Print() << "Writing " << stdDevOutfile << "..." << std::endl;
      WriteSingleLevelPlotfile(stdDevOutfile,*stdDevMF,names_stddev,geoms,timeFinal-time0,levelSteps);
      stdDevMF.reset();
    }
    
    }
    amrex::Finalize();
    return 0;
}



