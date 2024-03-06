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
    std::cerr << argv[0] << " infile=<name> outfile=<name> remove=<comps_name_to_remove>";
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


	// check if compiled as float
	{
	    Real test_float = 1.2;
	    if (test_float != float(test_float))
		amrex::Abort("WARNING: compiled as double");
	}
	
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
	    amrex::Abort("It would be safer not to have output name == input name");
	
	DataServices::SetBatchMode();
	Amrvis::FileType fileType(Amrvis::NEWPLT);

	// Set up for reading pltfile
	DataServices dataServices(infile, fileType);
	
	if (!dataServices.AmrDataOk())
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);
	
	AmrData& amrData = dataServices.AmrDataRef();

	
	// plot file variables
	const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
	const int init_ncomp = plotVarNames.size();

	// variables to remove provided by user
	const int remove_nComp = pp.countval("remove");
	Vector<std::string> removeVarNames;
	removeVarNames.resize(remove_nComp);
	pp.queryarr("remove",removeVarNames,0,remove_nComp);
	Vector<int> remove_comps;         // for index of remove variables
	
	// variables to keep	
	Vector<int> keep_comps;                // for index of keep variables
	Vector<std::string> keepVarNames; // for string of keep variables

	// initally set all variables to keep
	keepVarNames = plotVarNames;  
	for (int i = 0; i < plotVarNames.size(); ++i) {
	    keep_comps.push_back(i);
	}

	// find indexes to remove
	for (int i = 0; i < plotVarNames.size(); ++i) {
	    for (int j = 0; j < remove_nComp; ++j) {
		if (keepVarNames[i] == removeVarNames[j]) {
		    Print() << "removing plot variable: " << keepVarNames[i] << "\n";
		    remove_comps.push_back(i);
		}
	    }
	}

	// assign keep indexes to comps
	Vector<int> comps;
	for (const auto &value : keep_comps) {
	    if (std::find(remove_comps.begin(), remove_comps.end(), value) == remove_comps.end()) {
		comps.push_back(value);
	    }
	}	
	const int nComp = comps.size();
	
	
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
