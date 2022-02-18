#include <string>
#include <iostream>
//#include <set>
//#include <map>
#include <vector>

#include <AMReX_ParmParse.H>
//#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

typedef Array<Real,AMREX_SPACEDIM> dim3;
typedef Array<int,AMREX_SPACEDIM> int3;

#include <streamBinTubeStats.H>

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  if (ParallelDescriptor::NProcs()>1)
    Abort("Code is not yet parallel safe");
  
  ParmParse pp;

  // read infile name from inputs
  std::string infile;
  pp.get("infile",infile);

  // parse what files to write
  int writeStreamsToMatlab(0);
  pp.query("writeStreamsToMatlab",writeStreamsToMatlab);
  int writeTecplotSurfaceFromStream(0);
  pp.query("writeTecplotSurfaceFromStream",writeTecplotSurfaceFromStream);
  
  // declare size and data holders
  int nStreams, nElts, nPtsOnStream, nComps;
  std::vector<std::string> variableNames;
  Vector<int> faceData;
  Vector<Vector<Real>> streamData;

  // Read file
  readStreamBin(infile, nStreams, nElts, nPtsOnStream, nComps,
		variableNames, faceData, streamData);

  // report
  Print() << "nStreams      = " << nStreams << std::endl;
  Print() << "nElements     = " << nElts << std::endl;
  Print() << "nPtsOnStreams = " << nPtsOnStream << std::endl;
  Print() << "nComps        = " << nComps << std::endl;
  Print() << "Variable names:" << std::endl;
  for (int iComp=0; iComp<nComps; iComp++)
    Print() << "   " << iComp << ": " << variableNames[iComp] << std::endl;

  // let's write a matlab file for each variable
  if (writeStreamsToMatlab) {
    Print() << "Writing streams as matlab files..." << std::endl;
    writeStreamsMatlab(infile,nStreams,nPtsOnStream,nComps,variableNames,streamData);
  }
  
  // let's output a surface
  if (writeTecplotSurfaceFromStream) {
    Print() << "Writing a tecplot surface..." << std::endl;
    writeSurfaceFromStreamTecplot(infile,nStreams,nElts,nPtsOnStream,nComps,
				  variableNames,faceData,streamData);
  }

  // parse which variables to average
  int nAvg = pp.countval("avgComps");
  Vector<std::string> avgComps;
  Vector<int> avgIdx(nAvg);
  if (nAvg>0) {
    pp.getarr("avgComps",avgComps);
    Print() << "Average components:" << std::endl;
    for (int iAvg=0; iAvg<nAvg; iAvg++) {
      avgIdx[iAvg]=-1;
      for (int iComp=0; iComp<nComps; iComp++)
	if (variableNames[iComp]==avgComps[iAvg]) avgIdx[iAvg]=iComp;
      if (avgIdx[iAvg]==-1) {
	std::string msg="avgComp (" + avgComps[iAvg] + ") not in stream file";
	Abort(msg);
      } else {
	Print() << "   " << avgIdx[iAvg] << ": " << avgComps[iAvg] << std::endl;
      }
    }
  }

  // parse which variables to integrate
  int nInt = pp.countval("intComps");
  Vector<std::string> intComps;
  Vector<int> intIdx(nInt);
  if (nInt>0) {
    pp.getarr("intComps",intComps);
    Print() << "Average components:" << std::endl;
    for (int iInt=0; iInt<nInt; iInt++) {
      intIdx[iInt]=-1;
      for (int iComp=0; iComp<nComps; iComp++)
	if (variableNames[iComp]==intComps[iInt]) intIdx[iInt]=iComp;
      if (intIdx[iInt]==-1) {
	std::string msg="intComp (" + intComps[iInt] + ") not in stream file";
	Abort(msg);
      } else {
	Print() << "   " << intIdx[iInt] << ": " << intComps[iInt] << std::endl;
      }
    }
  }

  // Override for now
  nInt=0;
  
  // parse derived quantities (e.g. principal curvature zone)
  int nDer(0);
  Vector<std::string> derComps;
  
  // make space to hold output surface
  // for each element (i.e. triangle), we have three coordinates and one set of data
  // data will be written in triplicate, but no need to store all that
  // connectivity follows naturally from construction
  int nSurfComps(nAvg+nInt+nDer);
  Vector<Vector<dim3>> surfLocs(nElts);
  Vector<Vector<Real>> surfData(nElts);
  for (int iElt=0; iElt<nElts; iElt++) {
    surfLocs[iElt].resize(AMREX_SPACEDIM);
    surfData[iElt].resize(nSurfComps);    
  }

  // the surface is the mid point of the stream
  int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero
  
  // get the three stream indices from connectivity faceData
  Vector<int3> sIdx(nElts);
  for (int iElt=0; iElt<nElts; iElt++)
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      sIdx[iElt][iCorner] = faceData[iElt*AMREX_SPACEDIM+iCorner];
  
  // set locations
  Print() << "Setting locations ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    // define the spatial location of the three corners of the triangle
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) { // three corners (streams)
      int iStream=faceData[iElt*AMREX_SPACEDIM+iCorner];
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	surfLocs[iElt][iCorner][d] = streamData[iStream][nPtsOnStream*d+surfPt];
      }
    }
  }

  // evaluate volume
  Print() << "Evaluating volumes ..." << std::endl;
  Vector<Real> eltVol(nElts);
  for (int iElt=0; iElt<nElts; iElt++) {
    int3 iStream; // stream indices that make up this element
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      iStream[iCorner]=sIdx[iElt][iCorner];
    Real vol=0.;
    for (int iPt=1; iPt<nPtsOnStream; iPt++) {
      dim3 A, B, C, D, E, F;
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	A[d] = streamData[iStream[0]][nPtsOnStream*d+iPt-1];
	B[d] = streamData[iStream[1]][nPtsOnStream*d+iPt-1];
	C[d] = streamData[iStream[2]][nPtsOnStream*d+iPt-1];
	D[d] = streamData[iStream[0]][nPtsOnStream*d+iPt];
	E[d] = streamData[iStream[1]][nPtsOnStream*d+iPt];
	F[d] = streamData[iStream[2]][nPtsOnStream*d+iPt];
      }
      vol += wedge_volume(A,B,C,D,E,F);
    }
    eltVol[iElt]=vol;
  }

  // calculate averages
  Print() << "Calculating averages ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    // evaluate average of each component over the three points
    for (int iAvg=0; iAvg<nAvg; iAvg++) {
      int iComp=avgIdx[iAvg]; // component index that we're averaging
      Real avgVal(0.);
      for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	int iStream=sIdx[iElt][iCorner];
	avgVal += streamData[iStream][nPtsOnStream*iComp+surfPt];
      }
      surfData[iElt][iAvg] = avgVal/3.;
    }
  }

  // calculate integrals
  Print() << "Calculating integrals ..." << std::endl;

  // calculate derived
  Print() << "Calculating derived quanities ..." << std::endl;
  
  // write surface
  Print() << "Writing surface ..." << std::endl;
  writeSurfaceTecplot(infile, nAvg, avgComps, nInt, intComps, nDer, derComps,
		      nElts, eltVol, surfLocs, surfData);
  
}

//
// routine to read aja's binary stream files
//

void readStreamBin(std::string infile,
		   int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<int>& faceData,
		   Vector<Vector<Real>> &streamData)
{
  // Open header file
  std::string headerName = infile + "/Header";
  Print() << "Opening " << headerName << std::endl;
  std::ifstream ifs(headerName.c_str());
  std::istream* is = (infile=="-" ? (std::istream*)(&std::cin) : (std::istream*)(&ifs) );

  // read dummy header line
  std::string dummy;
  std::getline(ifs,dummy);

  // number of files to read
  int nFiles(-1);  
  ifs >> nFiles;
  Print() << "nFiles = " << nFiles << std::endl;

  // number of points on each stream
  ifs >> nStreams;
  Print() << "nStreams = " << nStreams << std::endl;

  // number of points on each stream
  ifs >> nPtsOnStream;
  Print() << "nPtsOnStream = " << nPtsOnStream << std::endl;

  // number of components
  ifs >> nComps;
  Print() << "nComps = " << nComps << std::endl;

  // next line
  std::getline(ifs,dummy);

  // read variable names
  variableNames.resize(nComps);
  variableNames = parseVarNames(*is);
  if (nComps!=variableNames.size())
    Abort("nComps != variableNames.size()");
  //Print() << "Variable names:" << std::endl;
  //for (int iComp=0; iComp<nComps; iComp++)
  //Print() << "   " << iComp << ": " << variableNames[iComp] << std::endl;

  // connectivity data
  int fds;
  ifs >> fds;
  Print() << "faceDataSize = " << fds << std::endl;
  std::getline(ifs,dummy);
  faceData.resize(fds);
  ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());

  nElts = fds/3;
  Print() << "nElts = " << nElts << std::endl;

  // close header
  ifs.close();

  //
  // give the streams a home
  //
  streamData.resize(nStreams+1);
  for (int iStream=0; iStream<=nStreams; iStream++)
    streamData[iStream].resize(nPtsOnStream*nComps);
  
  // keep fread happy
  size_t read_size;

  //
  // now loop over binary files
  //
  for (int iFile=0; iFile<nFiles; iFile++) {
    // read binary stream file
    std::string rootName = infile + "/str_";
    std::string fileName = Concatenate(rootName,iFile) + ".bin";
    //Print() << "Opening " << fileName << std::endl;
    FILE *file=fopen(fileName.c_str(),"r");

    int nFileStreams;
    read_size = fread(&(nFileStreams),sizeof(int),1,file);

    // loop over particle streams as written by pti (i.e. mangled order)
    for (int pindex=0; pindex<nFileStreams; pindex++) {
      // use the particle id to load data into right memory destination
      int iStream;
      read_size = fread(&(iStream),sizeof(int),1,file); // id
      
      for (int iComp=0; iComp<nComps; iComp++) {
	// by loading into [iStream] index, we're unmangling the parrallel particles
	int offset = iComp*nPtsOnStream;
	read_size = fread(&(streamData[iStream][offset]),sizeof(Real),nPtsOnStream,file);
      }      
    }
    
    fclose(file);

  } // iFiles

  Print() << "Finished reading stream data." << std::endl;
  
}

void
writeStreamsMatlab(std::string infile,
		   int& nStreams, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<Vector<Real>>& streamData)
{
  FILE *file;
  std::string filename;
  for (int iComp=0; iComp<nComps; iComp++) {
    filename = infile+"/"+variableNames[iComp]+".dat";
    file = fopen(filename.c_str(),"w");
    for (int iStream=1; iStream<=nStreams; iStream++) {
      for (int iPt=0; iPt<nPtsOnStream; iPt++) {
	fprintf(file,"%e ",streamData[iStream][nPtsOnStream*iComp+iPt]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  return;
}

void
writeSurfaceFromStreamTecplot(std::string infile,
			      int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
			      std::vector<std::string>& variableNames,
			      Vector<int>& faceData,
			      Vector<Vector<Real>>& streamData)
{
  std::string filename=infile+"_surfTec.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  std::string vars("VARIABLES =");
  for (int iComp=0; iComp<nComps; iComp++)
    vars += " " + variableNames[iComp];
  os << vars << std::endl;

  os << "ZONE T=\"streamBinTubeSurface\""
     << " N=" << nStreams
     << " E=" << nElts
     << " F=FEPOINT ET= TRIANGLE"
     << std::endl;

  int iPt=(nPtsOnStream-1)/2;
  for (int iStream=1; iStream<=nStreams; iStream++) {
    for (int iComp=0; iComp<nComps; iComp++) {
      os << streamData[iStream][nPtsOnStream*iComp+iPt] << " ";
    }
    os << std::endl;
  }

  for (int iElt=0; iElt<nElts; iElt++) {
    int offset=iElt*3;
    os << faceData[offset] << " "
       << faceData[offset+1] << " "
       << faceData[offset+2] << std::endl;
  }

  os.close();
  
  return;
}


void
writeSurfaceTecplot(std::string infile,
		    int& nAvg, Vector<std::string>& avgComps,
		    int& nInt, Vector<std::string>& intComps,
		    int& nDer, Vector<std::string>& derComps,
		    int& nElts, Vector<Real>& eltVol,
		    Vector<Vector<dim3>>& surfLocs,
		    Vector<Vector<Real>>& surfData)
{
  std::string filename=infile+"_binVolInt.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  std::string vars("VARIABLES = X Y Z volume");
  for (int iAvg=0; iAvg<nAvg; iAvg++)
    vars += " " + avgComps[iAvg] + "_avg";
  for (int iInt=0; iInt<nInt; iInt++)
    vars += " " + intComps[iInt] + "_volInt";
  for (int iDer=0; iDer<nDer; iDer++)
    vars += " " + derComps[iDer];
  os << vars << std::endl;

  os << "ZONE T=\"streamBinTubeSurface\""
     << " N=" << nElts*AMREX_SPACEDIM
     << " E=" << nElts
     << " F=FEPOINT ET=TRIANGLE"
     << std::endl;

  // write averages
  for (int iElt=0; iElt<nElts; iElt++) {
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
      // coordinate
      for (int d=0; d<AMREX_SPACEDIM; d++)
	os << surfLocs[iElt][iCorner][d] << " ";
      // volume
      os << eltVol[iElt] << " ";
      // averages
      for (int iAvg=0; iAvg<nAvg; iAvg++)
	os << surfData[iElt][iAvg] << " ";
      // integrals
      for (int iInt=0; iInt<nInt; iInt++)
	os << surfData[iElt][iInt] << " ";
      // derived
      for (int iDer=0; iDer<nDer; iDer++)
	os << surfData[iElt][iDer] << " ";
      os << std::endl;
    }
  }

  // write connectivity
  for (int iElt=0, idx=1; iElt<nElts; iElt++)
    os << idx++ << " " << idx++ << " " << idx++ << std::endl;

  os.close();
}

std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return amrex::Tokenize(line,std::string(" "));
}

Real
tetVol(const dim3& A, const dim3& B, const dim3& C, const dim3& D)
{
    Real R1[3], R2[3], R3[3], R4[3];

    for (int i=0; i<3; ++i)
    {
        R1[i] = D[i] - A[i];
        R2[i] = B[i] - A[i];
        R3[i] = C[i] - A[i];
    }

    for (int i=0; i<3; ++i)
    {
        R4[0] = R2[1]*R3[2] - R3[1]*R2[2];
        R4[1] = R2[2]*R3[0] - R3[2]*R2[0];
        R4[2] = R2[0]*R3[1] - R3[0]*R2[1];
    }
    
    Real res=0;
    for (int i=0; i<3; ++i)
        res += R1[i]*R4[i];

    return std::abs(res/6.);
}

Real
wedge_volume(const dim3& A, const dim3& B, const dim3& C,
	     const dim3& D, const dim3& E, const dim3& F)
{
    Real result;

    const Real vol_EABC = tetVol(A,B,C,E);
    const Real vol_ADEF = tetVol(A,D,E,F);
    const Real vol_ACEF = tetVol(C,E,F,A);
    
    result = (vol_EABC + vol_ADEF + vol_ACEF);

    return result;
}

Real
wedge_volume_int(const dim3& A, const Real vA,
		 const dim3& B, const Real vB,
		 const dim3& C, const Real vC,
		 const dim3& D, const Real vD,
		 const dim3& E, const Real vE,
		 const dim3& F, const Real vF)
{
    Real result;

    const Real vol_EABC = tetVol(A,B,C,E);
    const Real vol_ADEF = tetVol(A,D,E,F);
    const Real vol_ACEF = tetVol(C,E,F,A);
    const Real vol_DABC = tetVol(A,B,C,D);
    const Real vol_FABC = tetVol(A,B,C,F);
    const Real vol_BDEF = tetVol(B,D,E,F);
    const Real vol_CDEF = tetVol(C,D,E,F);
    const Real vol_ACED = tetVol(C,E,D,A);
    const Real vol_BCDF = tetVol(B,C,D,F);
    const Real vol_BCDE = tetVol(B,C,D,E);
    const Real vol_ABDF = tetVol(B,D,F,A);
    const Real vol_ABEF = tetVol(B,E,F,A);
    
    const Real int_1 = ( (vD+vA+vB+vC)*vol_DABC + 
			 (vB+vD+vE+vF)*vol_BDEF + 
			 (vB+vC+vD+vF)*vol_BCDF )/4.;
    
    const Real int_2 = ( (vD+vA+vB+vC)*vol_DABC + 
			 (vC+vD+vE+vF)*vol_CDEF + 
			 (vB+vC+vD+vE)*vol_BCDE )/4.;
    
    const Real int_3 = ( (vE+vA+vB+vC)*vol_EABC + 
			 (vA+vD+vE+vF)*vol_ADEF + 
			 (vA+vC+vE+vF)*vol_ACEF )/4.;
    
    const Real int_4 = ( (vE+vA+vB+vC)*vol_EABC + 
			 (vC+vD+vE+vF)*vol_CDEF + 
			 (vA+vC+vE+vD)*vol_ACED )/4.;
    
    const Real int_5 = ( (vF+vA+vB+vC)*vol_FABC + 
			 (vA+vD+vE+vF)*vol_ADEF + 
			 (vA+vB+vE+vF)*vol_ABEF )/4.;
    
    const Real int_6 = ( (vF+vA+vB+vC)*vol_FABC + 
			 (vB+vD+vE+vF)*vol_BDEF + 
			 (vA+vB+vD+vF)*vol_ABDF )/4.;
    
    result = (int_1 + int_2 + int_3 + int_4 + int_5 + int_6)/6.;

    return result;
}
