#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> [options] \n\tOptions:\n";
    std::cerr << "\t     outfile=<name>\n";
    exit(1);
}

static std::string parseTitle(std::istream& is);
static std::vector<std::string> parseVarNames(std::istream& is);
static std::string rootName(const std::string in);

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

    bool verbose = false; pp.query("verbose",verbose);
    
    std::string cVar;
    pp.get("cVar",cVar);
    Real cMin,cMax;
    pp.get("cMin",cMin);
    pp.get("cMax",cMax);
    std::string infile; pp.get("infile",infile);
    int cComp = -1; 
    
    std::ifstream ifs;
    std::istream* is = (infile=="-" ? (std::istream*)(&std::cin) : (std::istream*)(&ifs) );
    if (infile!="-")
    {
        if (verbose)
            std::cerr << "Opening " << infile << std::endl;
        ifs.open(infile.c_str(),std::ios::in|std::ios::binary);
        if (ifs.fail())
            std::cerr << "Unable to open file : " << infile << std::endl;
    }
    else
    {
        if (verbose)
            std::cerr << "Reading from stream" << std::endl;
    }
    const std::string title = parseTitle(*is);
    const std::vector<std::string> names = parseVarNames(*is);
    const int nComp = names.size();
    for (int i = 0; i<nComp; i++) {
      if (names[i] == cVar) {
	cComp = i;
	break;
      }
    }
    if (cComp < 0) {
      Abort("Conditional variable not in list");
    }
    int nElts;
    int MYLEN;
    (*is) >> nElts;
    (*is) >> MYLEN;

    FArrayBox nodeFab;
    nodeFab.readFrom((*is));
    Real* nodeData = nodeFab.dataPtr();
    int nPts = nodeFab.box().numPts();
    
    Vector<int> connData(nElts*MYLEN,0);
    (*is).read((char*)connData.dataPtr(),sizeof(int)*connData.size());

    // Build outfile name
    std::vector<std::string> infileTokens = amrex::Tokenize(infile,".");
    std::string outfile;
    outfile = infileTokens[0];
    for (int i=1; i<infileTokens.size()-1; ++i)
        outfile += std::string(".") + infileTokens[i];
    outfile += ".dat";
    pp.query("outfile",outfile);

    std::string vars("VARIABLES =");
    for (int j=0; j<names.size(); ++j)
        vars += " " + names[j];


    Vector<int> conditionedElt(nElts,0);
    Vector<int> keepPt(nPts,0);
    for (int i=0; nElts; ++i) {
      int offsetElt = i*MYLEN;
      int localTest = 0;
      for (int k = 0; k<MYLEN; k++) {
	int node = connData[offsetElt+k]-1;
	int nodeOffset = node*nComp;
	if (nodeData[nodeOffset+cComp] < cMax && nodeData[nodeOffset+cComp] > cMin) {
	  localTest += 1;
	}
      }
      if(localTest == 3) {
	conditionedElt[i] = 1;
	for (int k = 0; k<MYLEN; k++) {
	  int node = connData[offsetElt+k]-1;
	  keepPt[node] = 1;
	}
      }
    }
    int nCondPts = 0;
    int nCondElts = 0;
    for (int i = 0; i<nPts; i++) {
      nCondPts += keepPt[i];
    }
    for (int i = 0; i<nElts; i++) {
      nCondElts += conditionedElt[i];
    }

    Vector<int> ptMap(nPts,-1);
    int ptCount = 0;
    for (int i = 0; i < nPts; i++) {
      if(keepPt[i] == 1) {
	ptMap[ptCount] = i;
	ptCount++;
      }
    }
    // Write ASCII output file
    std::ofstream os(outfile.c_str(),std::ios::out);
    os << vars << std::endl;
    os << "ZONE T=\"" << title << "\" N=" << nCondPts << " E=" << nCondElts
       << " F=FEPOINT ET=";
    os << (MYLEN==2 ? "LINESEG" : "TRIANGLE") << std::endl;
    
    Vector<Real> condNodeData(nComp*nCondPts);
    Vector<int> condConnData(MYLEN*nCondElts);
    for (int i=0; i<nPts; ++i)
    {
      if (keepPt[i] == 1) {
        int offset = i*nComp;
        for (int k=0; k<nComp; ++k)
            os << nodeData[offset+k] << " ";
        os << std::endl;
      }
    }
    for (int i=0; i<nElts; ++i)
    {
      if (conditionedElt[i] == 1) {
	int offset = i*MYLEN;
        for (int k=0; k<MYLEN; ++k)
            os << ptMap[connData[offset+k]-1]+1 << " ";
        os << std::endl;
      }
    }
    os.close();
    
    }
    amrex::Finalize();
    return 0;
}

static
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return amrex::Tokenize(line,std::string(" "));
}

static std::string parseTitle(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return line;
}

static std::string rootName(const std::string inStr)
{
#ifdef WIN32
    const std::string dirSep("\\");
#else
    const std::string dirSep("/");
#endif
    std::vector<std::string> res = amrex::Tokenize(inStr,dirSep);
    res = amrex::Tokenize(res[res.size()-1],std::string("."));
    std::string out = res[0];
    for (int i=1; i<res.size()-1; ++i)
        out = out + std::string(".") + res[i];
    return out;
    //return (res.size()>1 ? res[res.size()-2] : res[res.size()-1]);
}

