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
/*
std::string
getFileRoot(const std::string& infile)
{
    Vector<std::string> tokens = Tokenize(infile,std::string("/"));
    return tokens[tokens.size()-1];
}
*/
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
    
    std::string outfile = "minPlot";
    pp.query("outfile",outfile);
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Print() << "Setting up minPlotfile..." << std::endl;
    // Set up for reading pltfile
    DataServices dataServices0(infiles[0], fileType);

    if (!dataServices0.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData0 = dataServices0.AmrDataRef();

    // read/write all comps   
    int nComp = pp.countval("comps");
    Vector<int> comps(nComp);
    if (nComp > 0) {
      pp.getarr("comps",comps);
      for (int i = 0; i<nComp; ++i) {
	Print() << "averaging comp " << comps[i] << std::endl;
      }
    } else {
      Print() << "averaging all comps" << std::endl;
      nComp = amrData0.NComp();
      comps.resize(nComp);
      for (int i=0; i<nComp; ++i) { 
	comps[i] = i;
      }
    }

    int finestLevel = amrData0.FinestLevel();
    //finest level to go down to
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Print() << "... averaging to resolution of level " << finestLevel << std::endl;
    // find probDomain at level we want
    Box probDomain0 = amrData0.ProbDomain()[finestLevel];
    //int max_grid_size = 32; pp.query("max_grid_size",max_grid_size);
    BoxArray ba=amrData0.boxArray(finestLevel);
    //tmpnewba.maxSize(max_grid_size);

    Print() << "... BoxArray set" << std::endl;
    Print() << "... number of boxes = " << ba.size() << std::endl;
    //const DistributionMapping& dm(newba);
    int levelSteps; //dummy variable
    
    RealBox rb(&(amrData0.ProbLo()[0]),
               &(amrData0.ProbHi()[0]));
    Vector<int> is_per(BL_SPACEDIM,0);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    int coord = 0;
    Geometry geoms(probDomain0,&rb,coord,&(is_per[0]));
    Print() << "... creating minMF" << std::endl;
    DistributionMapping dm(ba);
    MultiFab minMF(ba,dm,nComp,0);
    Print() << "... minMF created" << std::endl;
    Real time0 = amrData0.Time();
    Vector<std::string> names(nComp);
    Vector<int> destComps(nComp);
    for (int i=0; i<nComp; ++i) {
      names[i] = amrData0.PlotVarNames()[comps[i]];
      destComps[i]=i; 
      amrData0.FlushGrids(comps[i]); //do i need to flush if not filled?
    }
    Real timeFinal = 0;
    Print() << "minPlotfile setup completed." << std::endl;
    Long fab_megabytes = amrex::TotalBytesAllocatedInFabsHWM() / (1024*1024);
    ParallelDescriptor::ReduceLongMax(fab_megabytes, ParallelDescriptor::IOProcessorNumber());
    Print() << "Highest MFs mem. allocated on CPU (MB): " << fab_megabytes << std::endl;
    Print() << "Starting taking minimum... " << std::endl;
    for (int iFile = 0; iFile<numFiles; iFile++) {
      // Set up for reading pltfile
      DataServices dataServices(infiles[iFile], fileType);
      
      if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
      
      AmrData& amrData = dataServices.AmrDataRef();

      // make space
      MultiFab fileData(ba,dm,nComp,0); //should remove old multifab and replace with new one on loop
      Print() << "... loading " << infiles[iFile] << std::endl;
      // load data
      amrData.FillVar(fileData,finestLevel,names,destComps);
      Print() << "... " << infiles[iFile] << " loaded" << std::endl;
      Print() << "... flushing grids" << std::endl;
      for (int i = 0; i < nComp; i++) {
	  amrData.FlushGrids(comps[i]);
      }
      fab_megabytes = amrex::TotalBytesAllocatedInFabsHWM() / (1024*1024);
      ParallelDescriptor::ReduceLongMax(fab_megabytes, ParallelDescriptor::IOProcessorNumber());
      Print() << "... Highest MFs mem. allocated on CPU after loading fileData (MB): " << fab_megabytes << std::endl;
      Print() << "... finding local mins on MF" << std::endl;
      for (MFIter mfi(fileData,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  Array4<Real> const& inbox  = fileData.array(mfi);
	  Array4<Real> const& outbox = minMF.array(mfi);	      
	  AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
				  {
				    for (int n = 0; n < nComp; n++) {
				      if (iFile == 0) {
					outbox(i,j,k,n) = inbox(i,j,k,n);
				      } else {
					outbox(i,j,k,n) = std::min(inbox(i,j,k,n),outbox(i,j,k,n));
				      }							 
				    }
				    
				  });
	}
      if(iFile == numFiles-1) {
	timeFinal = amrData.Time();
      }
    }
    Print() << "Completed finding mins." << std::endl;
    Print() << "Writing " << outfile << "..." << std::endl;
    WriteSingleLevelPlotfile(outfile,minMF,names,geoms,timeFinal-time0,levelSteps);
    }
    amrex::Finalize();
    return 0;
}



