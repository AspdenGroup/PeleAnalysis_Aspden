#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AppendToPlotFile.H>

using namespace amrex;

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> outfile=<name> comps=<comp_names>";
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

	// get names of in ant outfiles
	std::string infile;
	pp.get("infile",infile);
	std::string outfile;
	pp.get("outfile",outfile);

	// check they are not the same (ie do not overwrite)
	if (infile == outfile)
	    amrex::Abort("it would be safer not to have output name == input name");
	
	DataServices::SetBatchMode();
	Amrvis::FileType fileType(Amrvis::NEWPLT);

	// Set up for reading pltfile
	DataServices dataServices(infile, fileType);
	
	if (!dataServices.AmrDataOk())
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);
	
	AmrData& amrData = dataServices.AmrDataRef();

	const int nComp = pp.countval("comps");
	Vector<std::string> scomps; scomps.resize(nComp);
	pp.getarr("comps",scomps,0,nComp);

	Vector<int> comps;
	const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
	for (int i = 0; i < plotVarNames.size(); ++i) {
	    for (int j = 0; j < scomps.size(); ++j) {
		if (plotVarNames[i] == scomps[j]) {
		    comps.push_back(i);
		}
	    }
	}
		
	int Nlev = amrData.FinestLevel() + 1;
	
	Vector<MultiFab*> fileData(Nlev);
	for (int lev=0; lev<Nlev; ++lev)
	{
	    const DistributionMapping dm(amrData.boxArray(lev));
	    fileData[lev] = new MultiFab(amrData.boxArray(lev),dm,nComp,0);
	}
	if (ParallelDescriptor::IOProcessor())
	    std::cerr << "Full MultiFab allocated " << std::endl;
	
	for (int lev=0; lev<Nlev; ++lev)
	{
	    for (int i=0; i<comps.size(); ++i)
	    {
		fileData[lev]->copy(amrData.GetGrids(lev,comps[i]),0,i,1);
		if (ParallelDescriptor::IOProcessor())
		    std::cerr << "After GetGrids: " << amrData.PlotVarNames()[comps[i]] << std::endl;
		if (ParallelDescriptor::IOProcessor())
		    std::cerr << "AmrData flushed: " << amrData.PlotVarNames()[comps[i]] << std::endl;
	    }
	}
	if (ParallelDescriptor::IOProcessor())
	    std::cerr << "File data loaded" << std::endl;
	
	Real progMin, progMax;
	for (int i=0; i<comps.size(); ++i) {
	    amrData.MinMax(amrData.ProbDomain()[Nlev-1], amrData.PlotVarNames()[comps[i]], Nlev-1, progMin, progMax);
	    if (ParallelDescriptor::IOProcessor()) {
		std::cout << amrData.PlotVarNames()[comps[i]] << " min/max: " << progMin << ", " << progMax << std::endl;
	    }
	}
	
	Vector<std::string> names(nComp);
	for (int i=0; i<comps.size(); ++i)
	    names[i] = amrData.PlotVarNames()[comps[i]];
//	for (int i=0; i<compsR.size(); ++i)
//	    names[compsL.size()+i] = amrDataR.PlotVarNames()[compsR[i]];
	
	bool verb = false;
	WritePlotFile(fileData,amrData,outfile,verb,names);
    }
    amrex::Finalize();
    return 0;
}
