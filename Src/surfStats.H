typedef Array<int, AMREX_SPACEDIM> int3;
typedef Array<Real,AMREX_SPACEDIM> dim3;

const Real half    = 1./2.;
const Real third   = 1./3.;
const Real quarter = 1./4.;
const Real sixth   = 1./6.;

std::vector<std::string> parseVarNames(std::istream& is);

int getVarIdx(std::string varName, Vector<std::string> variableNames, int nComps);


Real calcAvgVal(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData);

//Real calcMax(int compIdx, int nPtsOnStream, Vector<int3> sIdx, int iElt, Vector<Vector<Real>> streamData);

Real calcMax(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData);

//Real calcIntegral(int compIdx, int nPtsOnStream, Vector<int3> sIdx, int iElt, Vector<Vector<Real>> streamData, Vector<Real> eltArea);

Real calcIntegral(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData, Real eltArea);

void readStreamBin(std::string infile,
		   int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<int>& faceData,
		   Vector<Vector<Real>> &streamData);

void checkUnity1(Vector<Real> arr, int arrsize);

void checkUnity2(Vector<Vector<Real>> arr, int arrsize1, int arrsize2);

void writeSingleArray(Vector<Real> arr, int arrsize, std::string filename);

void writeDoubleArray(Vector<Vector<Real>> arr, int arrsize1, int arrsize2, std::string filename);

Real
wedge_area(const dim3& A,
	   const dim3& B,
	   const dim3& C);

Real
tetVol(const dim3& A, const dim3& B, const dim3& C, const dim3& D);

Real
wedge_volume(const dim3& A, const dim3& B, const dim3& C,
	     const dim3& D, const dim3& E, const dim3& F);

Real
wedge_volume_int(const dim3& A, const Real vA,
		 const dim3& B, const Real vB,
		 const dim3& C, const Real vC,
		 const dim3& D, const Real vD,
		 const dim3& E, const Real vE,
		 const dim3& F, const Real vF);
