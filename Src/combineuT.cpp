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

    std::string infileU; pp.get("infileU",infileU);
    std::string infileT; pp.get("infileT",infileT);
    
    std::string outfile = "TUcombine";
    pp.query("outfile",outfile);
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Print() << "Setting up plotfile..." << std::endl;
    // Set up for reading pltfile
    DataServices dataServicesu(infileU, fileType);
    DataServices dataServicest(infileT, fileType);
    
    if (!dataServicesu.AmrDataOk() || !dataServicest.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrDatau = dataServicesu.AmrDataRef();
    AmrData& amrDatat = dataServicest.AmrDataRef();

    int finestLevel = amrDatat.FinestLevel();
    //finest level to go down to
    pp.query("finestLevel",finestLevel);
    std::string uName = "z_velocity"; pp.query("uName",uName);
    std::string tName = "IDT"; pp.query("tName",tName);
    Box probDomain = amrDatat.ProbDomain()[finestLevel];
    BoxArray newba(probDomain);
    
    int levelSteps; //dummy variable
    
    RealBox rb(&(amrDatat.ProbLo()[0]),
               &(amrDatat.ProbHi()[0]));
    Vector<int> is_per(BL_SPACEDIM,0);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    int coord = 0;
    Geometry geoms(probDomain,&rb,coord,&(is_per[0]));
    std::unique_ptr<MultiFab> tuMF;
    int nComp = 3; //t, u, t*u
    tuMF.reset(new MultiFab(newba,DistributionMapping(newba),3,0));
    Vector<std::string> names(nComp);
    Vector<int> destComps(nComp);
    Vector<std::string> unames = {uName};
    Vector<std::string> tnames = {tName};
    Vector<std::string> outNames = {uName,tName,uName+tName};
    
    
    // load data
    Print() << "loading u" << std::endl;
    amrDatau.FillVar(*tuMF,finestLevel,unames,{0});
    Print() << "loading t" << std::endl;
    amrDatat.FillVar(*tuMF,finestLevel,tnames,{1});
    Print() << "multiplying" << std::endl;
    MultiFab::Copy(*tuMF,*tuMF,0,2,1,0);
    MultiFab::Multiply(*tuMF,*tuMF,1,2,1,0);
    Print() << "Writing " << outfile << "..." << std::endl;
    Real time = amrDatau.Time();
    WriteSingleLevelPlotfile(outfile,*tuMF,outNames,geoms,time,levelSteps);
    
    
    }
    amrex::Finalize();
    return 0;
}



