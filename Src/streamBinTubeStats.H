void readStreamBin(std::string infile,
		   int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<int>& faceData,
		   Vector<Vector<Real>> &streamData);

void
writeStreamsMatlab(std::string infile,
		   int& nStreams, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<Vector<Real>>& streamData);
  
void
writeSurfaceFromStreamTecplot(std::string infile,
			      int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
			      std::vector<std::string>& variableNames,
			      Vector<int>& faceData,
			      Vector<Vector<Real>>& streamData);
  
void
writeSurfaceTecplot(std::string infile,
		    int& nElts, Vector<Real>&        eltVol,   Vector<Vector<dim3>>& surfLocs,
		    int& nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		    int& nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		    int& nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer);

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

std::vector<std::string> parseVarNames(std::istream& is);

