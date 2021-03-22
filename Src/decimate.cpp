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

Real minLength(std::vector<std::vector<Real>> coords, std::vector<std::vector<int>> connData) {
  Real eltLength,x0,y0,x1,y1;
  int node0, node1;
  Real minimumLength = 1;
  for (int i = 0; i<connData[0].size(); i++) { //loop over elements
    node0 = connData[0][i]-1; //left node of element
    node1 = connData[1][i]-1; //right node of element
    x0 = coords[0][node0]; //x coord of left node
    y0 = coords[1][node0]; //y coord of left node
    x1 = coords[0][node1]; //x coord of right node
    y1 = coords[1][node1]; //y coord of right node
    eltLength = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)); //length of element
    if (eltLength < minimumLength) {
      minimumLength = eltLength;
    }
  }
  return minimumLength;
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

    bool verbose = false; pp.query("verbose",verbose);
    std::string infile;
    std::string outfile = "decimated.dat";
    Real minSize = 1;
    std::vector<double> prob_lo(2),prob_hi(2);
    std::vector<int> n_cell(2);
    int max_level;
    Real dx;
    pp.getarr("geometry.prob_lo",prob_lo,0,2);
    pp.getarr("geometry.prob_hi",prob_hi,0,2);
    pp.getarr("amr.n_cell",n_cell,0,2);
    pp.get("amr.max_level",max_level);
    //std::cout << "geometry= " << prob_lo[0] << " " << prob_lo[0]
    dx = (prob_hi[0]-prob_lo[0])/(double)(pow(2,max_level)*n_cell[0]);
    std::cout << "dx = " << dx << std::endl;
    
    pp.query("outfile",outfile);
    pp.get("infile",infile);
    pp.query("minSize",minSize);
    
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
    if (MYLEN > 2) {
      std::cerr << "use qslim" << std::endl;
    }
    FArrayBox nodeFab;
    nodeFab.readFrom((*is));
    Real* nodeData = nodeFab.dataPtr();
    int nPts = nodeFab.box().numPts();
    
    Vector<int> connData(nElts*MYLEN,0);
    (*is).read((char*)connData.dataPtr(),sizeof(int)*connData.size());

    std::vector<std::vector<Real>> coords(2, std::vector<Real>(nPts,0));
    std::vector<std::vector<int>> decConnData(2, std::vector<int>(nElts,0));
    for (int i = 0; i<nPts; i++) {
      coords[0][i] = nodeData[i*nComp];
      coords[1][i] = nodeData[i*nComp + 1];
    }
    for (int i = 0; i<nElts; i++) {
      decConnData[0][i] = connData[MYLEN*i];
      decConnData[1][i] = connData[MYLEN*i + 1];
    }
    int node0,node1,numDecElts,numDecPts;
    Real x0,y0,x1,y1,eltLength,xNew,yNew;

    std::cout << "minLength= " << minLength(coords,decConnData) << ", minSize= " << minSize*dx << std::endl;
    while(minLength(coords,decConnData) < minSize*dx) {
      for (int i = 0; i<decConnData[0].size(); i++) { //loop over elements
	numDecElts = decConnData[0].size(); //current number of elements
	numDecPts = coords[0].size();
	node0 = decConnData[0][i]-1; //left node of element
	node1 = decConnData[1][i]-1; //right node of element
	x0 = coords[0][node0]; //x coord of left node
	y0 = coords[1][node0]; //y coord of left node
	x1 = coords[0][node1]; //x coord of right node
	y1 = coords[1][node1]; //y coord of right node
	eltLength = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)); //length of element
	//std::cout << "eltLength: " << eltLength << ", minSize: " << minSize*dx << std::endl; 
	if (eltLength < minSize*dx) { //if too small
	  xNew = (x1+x0)/2; //create new x coord between the two
	  yNew = (y1+y0)/2; //create new y coord between the two
	  coords[0][node0] = xNew; //replace left node x coord with new x coord
	  coords[1][node0] = yNew; //replace right node y coord with new y coord
	  for (int j = i; j<numDecElts-1; j++) {
	    decConnData[0][j] = decConnData[0][j+1]; //remove left node from element list and shift list left
	    decConnData[1][j] = decConnData[1][j+1]; //remove right node from element list and shift list left
	  }
	  decConnData[0].resize(numDecElts-1); //resize - will remove duplicated last entry
	  decConnData[1].resize(numDecElts-1);
	  for (int j = 0; j<numDecElts-1; j++) {
	    if (decConnData[0][j] == node1+1) { //if left node of element is old node
	      decConnData[0][j] = node0+1; //attach element to new node
	    } else if (decConnData[1][j] == node1+1) {
	      decConnData[1][j] = node0+1; 
	    }
	    if (decConnData[0][j] > node1+1) {
	      decConnData[0][j] -= 1;
	    }
	    if (decConnData[1][j] > node1+1) {
	      decConnData[1][j] -= 1; //shift all right nodes down by one
	    }
	  }
	  for (int j = node1; j<nPts-1; j++) {
	    coords[0][j] = coords[0][j+1]; //remove right node x coord from list and shift left
	    coords[1][j] = coords[1][j+1]; //remove right node y coord from list and shift left
	  }
	  coords[0].resize(numDecPts-1); //resize - will remove duplicated last entry
	  coords[1].resize(numDecPts-1);
	  //break;
	}
      }
    }
    std::cout << "new minLength = " << minLength(coords,decConnData) << std::endl;

    numDecPts = coords[0].size();
    numDecElts = decConnData[0].size();
    std::string vars = "X Y";
    // Write ASCII output file
    std::ofstream os(outfile.c_str(),std::ios::out);
    os << "VARIABLES = ";
    os << vars << std::endl;
    os << "ZONE T=\"" << title << "\" N=" << numDecPts << " E=" << numDecElts
       << " F=FEPOINT ET=";
    os << (MYLEN==2 ? "LINESEG" : "TRIANGLE") << std::endl;
    for (int i = 0; i < numDecPts; i++) {
      os << coords[0][i] << " " << coords[1][i] << std::endl;
    }
    for (int i = 0; i < numDecElts; i++) {
      os << decConnData[0][i] << " " << decConnData[1][i] << std::endl;
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

