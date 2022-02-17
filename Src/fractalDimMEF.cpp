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
    std::string infile;
    std::string outfile;
    pp.get("infile",infile);
    
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
    std::string basename = infile.substr(0,18);	  
    outfile = basename+"_boxData.dat";
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
    
    Real maxX;
    Real minX;
    Real maxY;
    Real minY;
    Real maxXdiff;
    Real maxYdiff;
    Real minBoxSize;
    Real maxBoxSize;
    for (int i=0; i<nElts; ++i) {
      int offsetElt = i*MYLEN;
      if (MYLEN == 2) {
	int offsetData1 = connData[offsetElt]-1; //first node of element, -1 since node 1 corresponds to data at 0
	int offsetData2 = connData[offsetElt+1]-1; //second node of element
	Real x1 = nodeData[offsetData1*nComp]; //x coord of node 1
	Real x2 = nodeData[offsetData2*nComp]; //x coord of node 2
	Real y1 = nodeData[offsetData1*nComp+1]; //y coord of node 1
	Real y2 = nodeData[offsetData2*nComp+1]; //y coord of node 2
	if (i == 0) {
	  maxX = x1;
	  minX = x1;
	  maxY = y1;
	  minY = y1;
	}
	
	if (x1 > maxX) {
	  maxX = x1;
	}
	if (x1 < minX) {
	  minX = x1;
	}
	if (x2 > maxX) {
	  maxX = x2;
	}
	if (x2 < minX) {
	  minX = x2;
	}
	if (y1 > maxY) {
	  maxY = y1;
	}
	if (y1 < minY) {
	  minY = y1;
	}
	if (y2 > maxY) {
	  maxY = y2;
	}
	if (y2 < minY) {
	  minY = y2;
	}
	if ((maxX - minX) > (maxY-minY)) {
	  maxBoxSize = maxX - minX;
	} else {
	  maxBoxSize = maxY - minY;
	}
	
	
	if (sqrt((y2-y1)*(y2-y1)) > maxYdiff) {
	  maxYdiff = sqrt((y2-y1)*(y2-y1));
	}
	if (sqrt((x2-x1)*(x2-x1)) > maxXdiff) {
	  maxXdiff = sqrt((x2-x1)*(x2-x1));
	}
	if (maxXdiff > maxYdiff) {
	  minBoxSize = maxXdiff;
	} else {
	  minBoxSize = maxYdiff;
	}
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
	
	//NEED TO DO THIS FOR 3D ALSO
      }
    }
    pp.query("minBoxSize",minBoxSize);
    int subdivisions = int (maxBoxSize/minBoxSize);
    std::cout << "Subdivisions: " << subdivisions << std::endl;
    Real** dataArray = (Real**)malloc(subdivisions*sizeof(Real*));
    for (int s = 0; s<subdivisions; s++) {
      dataArray[s] = (Real*)malloc(2*sizeof(Real));
    }      
    //Real dataArray[subdivisions][2];
    for (int s = 1; s<=subdivisions; s++) {
      Real*** boxArray = (Real***)malloc(s*sizeof(Real**));
      for (int i = 0; i<s;i++) {
	boxArray[i] = (Real**)malloc(s*sizeof(Real*));
	for (int j = 0; j<s; j++) {
	  boxArray[i][j] = (Real*)malloc(3*sizeof(Real));
	}
      }      
      //Real boxArray[s][s][3];
      for (int i = 0; i<s; i++) {
	for (int j = 0; j<s; j++) {
	  boxArray[i][j][0] = minX + i*maxBoxSize/s; //x-coordinate
	  boxArray[i][j][1] = minY + j*maxBoxSize/s; //y-coordinate
	  boxArray[i][j][2] = 0; //number of elements in box
	}
      }
      for (int i = 0; i<nElts; i++) {
	int offsetElt = i*MYLEN;
	if (MYLEN == 2) {
	  int offsetData1 = connData[offsetElt]-1; //first node of element, -1 since node 1 corresponds to data at 0
	  int offsetData2 = connData[offsetElt+1]-1; //second node of element
	  Real x1 = nodeData[offsetData1*nComp]; //x coord of node 1
	  Real x2 = nodeData[offsetData2*nComp]; //x coord of node 2
	  Real y1 = nodeData[offsetData1*nComp+1]; //y coord of node 1
	  Real y2 = nodeData[offsetData2*nComp+1]; //y coord of node 2
	  int which_i_x1;
	  int which_i_x2;
	  int which_j_y1;
	  int which_j_y2;
	  if (x1 == maxX) {
	    which_i_x1 = s-1;
	  }
	  else {
	    which_i_x1 = int ((x1-minX)*s/maxBoxSize);
	  }
	  if (x2 == maxX) {
	    which_i_x2 = s-1;
	  } else {
	    which_i_x2 = int ((x2-minX)*s/maxBoxSize);
	  }
	  if (y1 == maxY) {
	    which_j_y1 = s-1;
	  } else {
	    which_j_y1 = int ((y1-minY)*s/maxBoxSize);
	  }
	  if (y2 == maxY) {
	    which_j_y2 = s-1;
	  } else { 
	    which_j_y2 = int ((y2-minY)*s/maxBoxSize);
	  }
	  if (which_i_x1 >= s || which_i_x2 >= s || which_j_y1 >= s || which_j_y2 >= s) { //debugging
	    std::cout << "Error: box index exceeds number of boxes" << std::endl;
	    std::cout << "number of boxes: " << s << std::endl;
	    std::cout << "node 1 box index: (" << which_i_x1 << "," << which_j_y1 << ")" << std::endl;
	    std::cout << "node 2 box index: (" << which_i_x2 << "," << which_j_y2 << ")" << std::endl;
	    exit(0);
	  }
	  if ((which_i_x1 == which_i_x2) && (which_j_y1 == which_j_y2)) { //element contained in one box
	    boxArray[which_i_x1][which_j_y1][2] += 1;
	  } else if ((which_i_x1 == which_i_x2) && (std::abs(which_j_y1-which_j_y2) == 1)) { //element crosses vertically
	    boxArray[which_i_x1][which_j_y1][2] += 1;
	    boxArray[which_i_x2][which_j_y2][2] += 1;
	  } else if ((which_j_y1 == which_j_y2) && (std::abs(which_i_x1-which_i_x2) == 1)) { //element cross horizontally
	    boxArray[which_i_x1][which_j_y1][2] += 1;
	    boxArray[which_i_x2][which_j_y2][2] += 1;
	  } else if ((std::abs(which_i_x1-which_i_x2) == 1) && (std::abs(which_j_y1-which_j_y2) == 1)) { //element crosses diagonally
	    boxArray[which_i_x1][which_j_y1][2] += 1;
	    boxArray[which_i_x2][which_j_y2][2] += 1;
	    Real grad = (y2-y1)/(x2-x1);
	    Real x_corner;
	    Real y_corner;
	    Real y_line;
	    if ((which_i_x1 == which_i_x2 + 1) && grad < 0) { //upper left quadrant
	      x_corner = boxArray[which_i_x1][which_j_y1][0];
	      y_corner = boxArray[which_i_x1][which_j_y1][1];
	      y_line = grad*(x_corner-x1) + y1;
	      if (y_line < y_corner) {
		boxArray[which_i_x2][which_j_y1][2] += 1;
	      } else {
		boxArray[which_i_x1][which_j_y2][2] += 1;
	      }
	    } else if ((which_i_x1 == which_i_x2 + 1) && grad > 0) { //lower left quadrant
	      x_corner = boxArray[which_i_x1][which_j_y2][0];
	      y_corner = boxArray[which_i_x1][which_j_y2][1];
	      y_line = grad*(x_corner-x1) + y1;
	      if (y_line < y_corner) {
		boxArray[which_i_x1][which_j_y2][2] += 1;
	      } else {
		boxArray[which_i_x2][which_j_y1][2] += 1;
	      }
	    } else if ((which_i_x1 == which_i_x2 - 1) && grad > 0) { //upper right quadrant
	      x_corner = boxArray[which_i_x2][which_j_y1][0];
	      y_corner = boxArray[which_i_x2][which_j_y1][1];
	      y_line = grad*(x_corner-x1) + y1;
	      if (y_line < y_corner) {
		boxArray[which_i_x2][which_j_y1][2] += 1;
	      } else {
		boxArray[which_i_x1][which_j_y2][2] += 1;
	      } 
	    } else if ((which_i_x1 == which_i_x2 -1) && grad < 0) { //lower right quadrant
	      x_corner = boxArray[which_i_x2][which_j_y2][0];
	      y_corner = boxArray[which_i_x2][which_j_y2][1];
	      y_line = grad*(x_corner-x1) + y1;
	      if (y_line < y_corner) {
		boxArray[which_i_x1][which_j_y2][2] += 1;
	      } else {
		boxArray[which_i_x2][which_j_y1][2] += 1;
	      }
	    } else {
	      std::cout << "dodgy diagonal" << std::endl;
	      exit(0);
	    }
	  } 
	} else {
	  std::cout << "can't detect what type of element box crossing" << std::endl;
	  exit(0);
	}	    	    
      }
      int nBoxes = 0;
      for (int i=0; i<s;i++) {
	for (int j=0; j<s;j++) {
	  if(boxArray[i][j][2] > 0) {
	    nBoxes += 1;
	  }
	}
      }
      free(boxArray);
      dataArray[s-1][0] = maxBoxSize/s;
      dataArray[s-1][1] = nBoxes;
    }
    
    
    //write ASCII file
    std::ofstream file(outfile.c_str(),std::ios::out);
    for (int s = 0; s<subdivisions;s++) {
      file << dataArray[s][0] << " " << dataArray[s][1] << std::endl;
    }
    file.close();
    free(dataArray);
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

