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
    
    std::string outfile = "flucPlot";
    pp.query("outfile",outfile);
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Print() << "Setting up avgPlotfile..." << std::endl;
    // Set up for reading pltfile
    DataServices dataServices0(infiles[0], fileType);

    if (!dataServices0.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData0 = dataServices0.AmrDataRef();
    
    // read/write comps

    //lets hard code for now
    // we want rho (1), rhou_i (2), rhou_iu_j (3), rhoY (3), rhoYu_i (6), rhoh (1), rhohu_i (2) 
    int nCompIn = 1+3+3+1;//pp.countval("vars");
    int nCompAvg = 1+2+3+3+6+1+2;
    //ui'uj' (3) Y'u' (6) h'u' (2)
    int nCompFluc = 3+6+2;
    int nCompOut = nCompAvg+nCompFluc;
    Vector<std::string> vars = {"density","x_velocity","y_velocity","z_velocity","Y(H2)","Y(O2)","Y(N2)","rhoh"};
    Vector<int> destfillcomps(nCompIn);
    Print() << "Averaging : " << std::endl;
    for (int i = 0; i<nCompIn; i++) {
      //pp.get("vars",vars[i],i);
      destfillcomps[i] = i;
      Print() << vars[i] << std::endl;
    }
    
    const Vector<std::string>& plotVarNames = amrData0.PlotVarNames();

    Vector<int> comps(nCompIn);
    for (int i=0; i<nCompIn; i++) {
      comps[i] = -1;
      for (int j=0; j<plotVarNames.size(); j++) {
	if (plotVarNames[j] == vars[i]) {
	  comps[i] = j;
	  break;
	}
      }
    }
    for (int i=0; i<nCompIn; i++) {
      if (comps[i] < 0) {
	Print() << vars[i] << " not found" << std::endl;
      }
    }


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
    avgMF.reset(new MultiFab(newba,DistributionMapping(newba),nCompOut,0));
    avgMF->setVal(0);
    Print() << "... avg MF created" << std::endl;
    Real time0 = amrData0.Time();
    Vector<std::string> names_avg = {"rho","ur","uz","YH2","YO2","YN2","h","ur2","uruz","uz2","YH2uz","YO2uz","YN2uz","huz","YH2ur","YO2ur","YN2ur","hur"};
    for (int i=0; i<nCompIn; ++i) {
      amrData0.FlushGrids(comps[i]); //do i need to flush if not filled?
    }
    Real previousTime = time0;
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
      fileData.reset(new MultiFab(newba,DistributionMapping(newba),nCompIn,0));
      Print() << "... loading " << infiles[iFile] << std::endl;
      // load data
      amrData.FillVar(*fileData,finestLevel,amrData.PlotVarNames(),comps);
      for (int i = 0; i < nCompIn; i++) {
	  amrData.FlushGrids(comps[i]);
      }
      
      std::unique_ptr<MultiFab> tmpMF;
      tmpMF.reset(new MultiFab(newba,DistributionMapping(newba),nCompAvg,0));
      for (MFIter mfi(*tmpMF,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	Array4<Real> fileArray = fileData->array(mfi);
	Array4<Real> tmpArray = tmpMF->array(mfi);
	AMREX_PARALLEL_FOR_3D (bx, i, j, k,
			       {
				 tmpArray(i,j,k,0) = fileArray(i,j,k,0); //density
				 tmpArray(i,j,k,1) = fileArray(i,j,k,0)*std::sqrt(fileArray(i,j,k,1)*fileArray(i,j,k,1) + fileArray(i,j,k,2)*fileArray(i,j,k,2)); //rhour
				 for (int n = 0; n<4; n++) {
				   tmpArray(i,j,k,2+n) = fileArray(i,j,k,0)*fileArray(i,j,k,3+n); //rhouz rhoYs
				 }
				 tmpArray(i,j,k,6) = fileArray(i,j,k,7); //rhoh
				 tmpArray(i,j,k,7) = tmpArray(i,j,k,1)*tmpArray(i,j,k,1)/tmpArray(i,j,k,0); //rhour^2
				 tmpArray(i,j,k,8) = tmpArray(i,j,k,1)*fileArray(i,j,k,3); //rhouruz
				 for (int n=0;n<4;n++) {
				   tmpArray(i,j,k,9+n) = tmpArray(i,j,k,2)*fileArray(i,j,k,3+n); //rhouz^2, rhoYuz
				 }
				 tmpArray(i,j,k,13) = fileArray(i,j,k,3)*fileArray(i,j,k,7); //rhohuz
				 for (int n=0;n<3;n++) {
				   tmpArray(i,j,k,14+n) = tmpArray(i,j,k,1)*fileArray(i,j,k,4+n); //rhoYur
				 }
				 tmpArray(i,j,k,17) = tmpArray(i,j,k,1)*fileArray(i,j,k,7)/fileArray(i,j,k,0); //rhohur				 
			       });
      }
      
      Long fab_megabytes = amrex::TotalBytesAllocatedInFabsHWM() / (1024*1024);
      ParallelDescriptor::ReduceLongMax(fab_megabytes, ParallelDescriptor::IOProcessorNumber());
      Print() << "... Highest MFs mem. allocated on CPU after loading fileData (MB): " << fab_megabytes << std::endl;

      Print() << "... " << infiles[iFile] << " loaded" << std::endl;
      Print() << "... adding to averaged MF" << std::endl;
      MultiFab::Saxpy(*avgMF,1.0/numFiles,*tmpMF,0,0,nCompAvg,0);
      if(iFile == numFiles-1) {
	timeFinal = amrData.Time();
      }      
    }
    for (int n = 1; n<nCompOut;n++) {
      MultiFab::Divide(*avgMF,*avgMF,0,n,1,0); //favre average
    }
    Print() << "Completed averaging." << std::endl;
    Print() << "Doing fluctuations" << std::endl;
    Print() << "... creating flucMF" << std::endl;
   
    std::unique_ptr<MultiFab> flucMF;
    flucMF.reset(new MultiFab(newba,DistributionMapping(newba),nCompFluc,0));
    flucMF->setVal(0);
    Print() << "... flucMF created" << std::endl;
    Vector<std::string> names_fluc = {"rho uz'uz'","rho uz'ur'","rho ur'ur'","rho YH2'uz'","rho YO2'uz'","rho YN2'uz'","rhoh'uz'","rho YH2'ur'","rho YO2'ur'","rho YN2'ur'","rho h'ur'"};
    Print() << " setup completed." << std::endl;
    
    for (MFIter mfi(*avgMF,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> avgArray = avgMF->array(mfi,0);
      Array4<Real> flucArray = avgMF->array(mfi,nCompAvg);
      AMREX_PARALLEL_FOR_3D (bx, i, j, k,
			     {
			       flucArray(i,j,k,0) = avgArray(i,j,k,0)*(avgArray(i,j,k,9)-avgArray(i,j,k,2)*avgArray(i,j,k,2)); //rhouzuz
			       flucArray(i,j,k,1) = avgArray(i,j,k,0)*(avgArray(i,j,k,8)-avgArray(i,j,k,2)*avgArray(i,j,k,1)); //rhouzur
			       flucArray(i,j,k,2) = avgArray(i,k,k,0)*(avgArray(i,j,k,7)-avgArray(i,j,k,1)*avgArray(i,j,k,1)); //rhourur
			       flucArray(i,j,k,3) = avgArray(i,j,k,0)*(avgArray(i,j,k,10)-avgArray(i,j,k,2)*avgArray(i,j,k,3)); //rhouzYH2
			       flucArray(i,j,k,4) = avgArray(i,j,k,0)*(avgArray(i,j,k,11)-avgArray(i,j,k,2)*avgArray(i,j,k,4)); //rhouzYO2
			       flucArray(i,j,k,5) = avgArray(i,j,k,0)*(avgArray(i,j,k,12)-avgArray(i,j,k,2)*avgArray(i,j,k,5)); //rhouzYN2
			       flucArray(i,j,k,6) = avgArray(i,j,k,0)*(avgArray(i,j,k,13)-avgArray(i,j,k,2)*avgArray(i,j,k,6)); //rhouzh
			       flucArray(i,j,k,7) = avgArray(i,j,k,0)*(avgArray(i,j,k,14)-avgArray(i,j,k,1)*avgArray(i,j,k,3)); //rhourYH2
			       flucArray(i,j,k,8) = avgArray(i,j,k,0)*(avgArray(i,j,k,15)-avgArray(i,j,k,1)*avgArray(i,j,k,4)); //rhourYO2
			       flucArray(i,j,k,9) = avgArray(i,j,k,0)*(avgArray(i,j,k,16)-avgArray(i,j,k,1)*avgArray(i,j,k,5)); //rhourYN2
			       flucArray(i,j,k,10) = avgArray(i,j,k,0)*(avgArray(i,j,k,17)-avgArray(i,j,k,1)*avgArray(i,j,k,6)); //rhourh	 
			     });
    }
    Vector<std::string> names(nCompOut);
    for (int n=0;n<nCompAvg;n++) {
      names[n] = names_avg[n];
    }
    for (int n=nCompAvg;n<nCompOut;n++) {
      names[n] = names_fluc[n-nCompAvg];
    }
    Print() << "Finished calculating flucs." << std::endl;
    WriteSingleLevelPlotfile(outfile,*avgMF,names,geoms,timeFinal-time0,levelSteps);
    }
    amrex::Finalize();
    return 0;
}



