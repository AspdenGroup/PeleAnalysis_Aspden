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
    int nSurfFiles(pp.countval("infile"));
    std::string infile;
    std::string outfile = "area.dat";
    pp.query("outfile",outfile);
    Real dataArr[nSurfFiles];
    std::string timeArr[nSurfFiles];
    for(int iSurf = 0; iSurf < nSurfFiles; ++iSurf) {
      pp.get("infile",infile,iSurf);
    
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

      Real totLength = 0;
      Real totArea = 0;
      for (int i=0; i<nElts; ++i)
	{
	  int offsetElt = i*MYLEN;
	  if (MYLEN == 2) {
	    int offsetData1 = connData[offsetElt]-1; //first node of element, -1 since node 1 corresponds to data at 0
	    int offsetData2 = connData[offsetElt+1]-1; //second node of element
	    Real x1 = nodeData[offsetData1*nComp]; //x coord of node 1
	    Real x2 = nodeData[offsetData2*nComp]; //x coord of node 2
	    Real y1 = nodeData[offsetData1*nComp+1]; //y coord of node 1
	    Real y2 = nodeData[offsetData2*nComp+1]; //y coord of node 2
	    //std::cout << "(" << x1 << "," << y1 << "),(" << x2 << "," << y2 << ")" << std::endl; //prints all coords (check)
	    Real eltLength = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	    totLength += eltLength;
	  } else {
	    int offsetData1 = connData[offsetElt]-1; //first node of element, -1 since node 1 corresponds to data at 0
	    int offsetData2 = connData[offsetElt+1]-1; //second node of element
	    int offsetData3 = connData[offsetElt+2]-1; //third node of element (since triangular)
	    Real x1 = nodeData[offsetData1*nComp]; //x coord of node 1
	    Real x2 = nodeData[offsetData2*nComp]; //x coord of node 2
	    Real x3 = nodeData[offsetData3*nComp]; //x coord of node 3
	    Real y1 = nodeData[offsetData1*nComp+1]; //y coord of node 1
	    Real y2 = nodeData[offsetData2*nComp+1]; //y coord of node 2
	    Real y3 = nodeData[offsetData3*nComp+1]; //y coord of node 3
	    Real z1 = nodeData[offsetData1*nComp+2]; //z coord of node 1
	    Real z2 = nodeData[offsetData2*nComp+2]; //z coord of node 2
	    Real z3 = nodeData[offsetData3*nComp+2]; //z coord of node 3
	    Real edgeLength1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)); //length of edge 1 -> 2
	    Real edgeLength2 = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2)); //length of edge 2 -> 3
	    Real edgeLength3 = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3)); //length of edge 3 -> 1
	    Real s = (edgeLength1+edgeLength2+edgeLength3)/2; //semi-perimeter of triangle
	    Real eltArea = sqrt(s*(s-edgeLength1)*(s-edgeLength2)*(s-edgeLength3)); //Heron's formula for area
	    totArea += eltArea;
	  }
	}
      if (MYLEN == 2) {
	dataArr[iSurf] = totLength;
      } else {
	dataArr[iSurf] = totArea;
      }
      timeArr[iSurf] = title;
    }
    //Write ASCII file
    std::ofstream os(outfile.c_str(),std::ios::out);
    for (int iSurf = 0; iSurf<nSurfFiles; ++iSurf) {
      os << timeArr[iSurf] << " " << dataArr[iSurf] << std::endl;
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

