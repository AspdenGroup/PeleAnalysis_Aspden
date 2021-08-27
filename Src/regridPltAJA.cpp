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

    std::string infile;
    pp.get("infile",infile);

    std::string outfile;
    pp.get("outfile",outfile);
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Set up for reading pltfile
    DataServices dataServices(infile, fileType);

    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData = dataServices.AmrDataRef();

    // query number of boxes to split finest level
    // defaul the number of cores
    int nBoxes(ParallelDescriptor::NProcs());
    pp.query("nBoxes",nBoxes);

    // read/write all comps
    int nComp = amrData.NComp();
    Vector<int> comps(nComp);
    for (int i=0; i<nComp; ++i)
      comps[i] = i;

    int finestLevel = amrData.FinestLevel();
    int Nlev = finestLevel + 1;
    pp.query("finestLevel",finestLevel);
    Nlev = std::max(0, std::min(Nlev, finestLevel+1));

    Vector<MultiFab*> fileData(Nlev);
    Vector<BoxArray>  newba(Nlev);
    Vector<DistributionMapping> dm(Nlev);
    
    for (int lev=0; lev<Nlev; ++lev) {
      // file box array
      BoxArray fileba  = amrData.boxArray(lev);

      if (lev==finestLevel) {
	// create a single box that covers the finest level
	// and divide it into the number of boxes specified
	Box      minBox = fileba.minimalBox();
	BoxList  newbl(minBox,nBoxes);
	BoxArray balev(newbl);
	DistributionMapping dmlev(balev);
	newba[lev]=balev;
	dm[lev]=dmlev;
	// some analytics
	if (ParallelDescriptor::IOProcessor()) {
	  int fileCount(0);
	  int nGrow(50);
	  pp.query("nGrow",nGrow);
	  for (int iBox=0; iBox<fileba.size(); iBox++) {
	    Box box=fileba[iBox];
	    fileCount += (box.size()[0]+nGrow)+(box.size()[1]+nGrow)+(box.size()[2]+nGrow);
	  }
	  int newCount(0);
	  for (int iBox=0; iBox<balev.size(); iBox++) {
	    Box box=balev[iBox];
	    newCount += (box.size()[0]+nGrow)+(box.size()[1]+nGrow)+(box.size()[2]+nGrow);
	  }
	  Real saving = ( (Real)newCount ) / ( (Real)fileCount );
	  std::cout << "Estimated saving = "
		    << newCount << " / " << fileCount << " = "
		    << saving << std::endl;
	}
      } else {
	// just use file boxArray for lower levels
	DistributionMapping dmlev(fileba);
	newba[lev]=fileba;
	dm[lev]=dmlev;
      }

      // make space
      fileData[lev] = new MultiFab(newba[lev],dm[lev],nComp,0);
      
      // load data
      amrData.FillVar(*fileData[lev],lev,amrData.PlotVarNames(),comps);
    }

    // flushGrids
    for (int lev=0; lev<Nlev; ++lev) 
        for (int i=0; i<comps.size(); ++i) 
            amrData.FlushGrids(comps[i]);

    Vector<std::string> names(nComp);
    for (int i=0; i<comps.size(); ++i)
        names[i] = amrData.PlotVarNames()[comps[i]];

    Vector<Geometry> geoms(Nlev);
    Vector<int> levelSteps(Nlev);
    Vector<IntVect> refRatio(Nlev-1);
    Vector<const MultiFab*> dat(Nlev);

    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    int coord = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      levelSteps[lev] = 666;
      if (lev < Nlev-1) {
        int r = amrData.RefRatio()[lev];
        refRatio[lev] = IntVect(D_DECL(r,r,r));
      }
      dat[lev] = fileData[lev];
    }
    
    WriteMultiLevelPlotfile(outfile,Nlev,dat,names,geoms,amrData.Time(),levelSteps,refRatio);
    }
    amrex::Finalize();
    return 0;
}






































#if 0
    //
    // AJA: start
    //
    int nAcross, nLayers;
    pp.get("nAcross",nAcross);
    pp.get("nLayers",nLayers);
    int nBoxes = nAcross*nAcross*nLayers;

    int nProcs  = ParallelDescriptor::NProcs();
    int myProcs = ParallelDescriptor::MyProc();

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Regridding to:" << '\n'
		<< "   nAcross = " << nAcross << '\n'
		<< "   nLayers = " << nLayers << '\n'
		<< "   hence nBoxes/nCores = "
		<< nBoxes << " / "
		<< nProcs << std::endl;
    }
    
    if (nBoxes!=nProcs)
      Abort("box-proc mismatch");
    
    //
    // AJA: end
    //
#endif




    //
    // AJA: start
    //
#if 0
    // Get a minimal box at the finest level to slice and dice
    int fLev = finestLevel;
   
    // take i and j from domain
    Box probDomain = amrData.ProbDomain()[fLev];

    int iMin = probDomain.smallEnd()[0];
    int iMax = probDomain.bigEnd()[0];
    int iLen = probDomain.size()[0];
    
    int jMin = probDomain.smallEnd()[1];
    int jMax = probDomain.bigEnd()[1];
    int jLen = probDomain.size()[1];

    // find k that covers the whole finest level
    BoxArray ba     = amrData.boxArray(fLev);
    Box      minBox = ba.minimalBox();

    int kMin = minBox.smallEnd()[2];
    int kMax = minBox.bigEnd()[2];
    int kLen = minBox.size()[2];
    
    if (!iLen%nAcross) Abort("Doesn't divide in x");
    if (!jLen%nAcross) Abort("Doesn't divide in y");

    int iGridSize = iLen/nAcross;
    int jGridSize = jLen/nAcross;
    int kGridSize = kLen/nLayers;
    int kRem      = kLen%nLayers;

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "probDomain = " << probDomain << std::endl;
      std::cout << "minBox = " << minBox << std::endl;
      std::cout << "i Min/Max/Len/GridSize = "
		<< iMin << " / "
		<< iMax << " / "
		<< iLen << " / "
		<< iGridSize << std::endl;
      std::cout << "j Min/Max/Len = "
		<< jMin << " / "
		<< jMax << " / "
		<< jLen << " / "
		<< jGridSize << std::endl;
      std::cout << "k Min/Max/Len = "
		<< kMin << " / "
		<< kMax << " / "
		<< kLen << " / "
		<< kGridSize << " ( "
		<< kRem << " ) " << std::endl;
    }
    
    int iBox=0;
    for (int k=0, bak=kMin; k<nLayers; k++) {

      int bakLen=kGridSize;
      if (k<kRem) bakLen++;

      int bakLo=bak;
      int bakHi=bak+bakLen-1;

      if (ParallelDescriptor::IOProcessor())
	std::cout << "Layer / lo / hi = "
		  << k << " / "
		  << bakLo << " / "
		  << bakHi << std::endl;

      for (int j=0; j<nAcross; j++) {
	int bajLo =  j    * jGridSize;
	int bajHi = (j+1) * jGridSize - 1;

	for (int i=0; i<nAcross; i++) {
	  int baiLo =  i    * iGridSize;
	  int baiHi = (i+1) * iGridSize - 1;

	  if (ParallelDescriptor::IOProcessor())
	    std::cout << "Box " << iBox << " : " 
		      << "(" << baiLo << "," << bajLo << "," << bakLo << "), "
		      << "(" << baiHi << "," << bajHi << "," << bakHi << ")"
		      << std::endl;
	  
	  iBox++;
	}
      }
      bak += bakLen;
    }

    BoxList boxList(minBox,nBoxes);
    std::cout << boxList << std::endl;
#endif
    //
    // AJA: end
    //
