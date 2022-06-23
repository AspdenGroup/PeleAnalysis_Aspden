#include <string>
#include <iostream>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

#include <streamBinTubeStats_TLH.H>

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
  int writeBasic(0);
  pp.query("writeBasic",writeBasic);
  std::string fuelName="H2";
  pp.query("fuelName",fuelName);
  // declare size and data holders
  int nStreams, nElts, nPtsOnStream, nComps;
  Vector<std::string> variableNames;
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
      avgIdx[iAvg] = getVarIdx(avgComps[iAvg],variableNames,nComps);
      Print() << "   " << avgIdx[iAvg] << ": " << avgComps[iAvg] << std::endl;
    }
  }

  // parse which variables to integrate
  int nInt = pp.countval("intComps");
  Vector<std::string> intComps;
  Vector<int> intIdx(nInt);
  if (nInt>0) {
    pp.getarr("intComps",intComps);
    Print() << "Integral components:" << std::endl;
    for (int iInt=0; iInt<nInt; iInt++) {
      intIdx[iInt] = getVarIdx(intComps[iInt],variableNames,nComps);
      Print() << "   " << intIdx[iInt] << ": " << intComps[iInt] << std::endl;
    }
  }

  int localDerOut = -1;
  // parse derived quantities (e.g. principal curvature zone)
  int nDerFlag = pp.countval("derComps");
  Vector<std::string> derCompsIn;
  Vector<std::string> derCompsOut;
  Vector<Vector<int>> derIdxIn(nDerFlag); //derIn
  Vector<Vector<int>> derIdxOut(nDerFlag); //derOut
  int nDerOut = 0;
  Real deltaT,rhoY,pkzFF; 
  if (nDerFlag>0) pp.getarr("derComps",derCompsIn);
  for (int iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
    if (derCompsIn[iDerFlag]=="flameThickness") {
      Real reacTemp, prodTemp;
      pp.get("reacTemp",reacTemp);
      pp.get("prodTemp",prodTemp);
      deltaT=prodTemp-reacTemp;
      std::string tempGradVar="ModGradTemp";
      pp.query("tempGradVar",tempGradVar);
      derIdxIn[iDerFlag].push_back(getVarIdx(tempGradVar,variableNames,nComps)); //add gradT idx
      derIdxOut[iDerFlag].resize(1); //one out idx (thermalthickness)
      localDerOut += 1;
      derIdxOut[iDerFlag][0] = localDerOut;
      derCompsOut.push_back("thermalThickness");
      nDerOut += 1;
    }
    if (derCompsIn[iDerFlag]=="flameSpeed") {
      pp.get("rhoY",rhoY);
      std::string FCRVar=fuelName+"_ConsumptionRate";
      pp.query("FCRVar",FCRVar);
      derIdxIn[iDerFlag].push_back(getVarIdx(FCRVar,variableNames,nComps)); //add FCR idx
      derIdxOut[iDerFlag].resize(1); //one out idx (flamespeed)
      localDerOut += 1;
      derIdxOut[iDerFlag][0] = localDerOut;
      derCompsOut.push_back("flameSpeed");
      nDerOut += 1;
    }
    if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
      std::string pkzMkVar = "MeanCurvature_prog_"+fuelName;
      std::string pkzGkVar = "GaussianCurvature_prog_"+fuelName;
      pp.query("pkzMkVar",pkzMkVar);
      pp.query("pkzGkVar",pkzGkVar);
      Real pkzLength;
      pp.get("pkzLength",pkzLength);
      pkzFF = 1.0 / (2*pkzLength);
      derIdxIn[iDerFlag].push_back(getVarIdx(pkzMkVar,variableNames,nComps)); //add mk idx
      derIdxIn[iDerFlag].push_back(getVarIdx(pkzGkVar,variableNames,nComps)); //add gk idx
      derIdxOut[iDerFlag].resize(3); //three out idx (zone,k1,k2)
      for (int iDer = 0; iDer<3; iDer++) {
	localDerOut += 1;
	derIdxOut[iDerFlag][iDer] = localDerOut;
      }
      derCompsOut.push_back("k1");
      derCompsOut.push_back("k2");
      derCompsOut.push_back("zone");
      nDerOut += 3;      
    }
  }
  AMREX_ALWAYS_ASSERT(derCompsOut.size() == nDerOut); //just a check on derive components out
  
  
  // make space to hold output surface
  // for each element (i.e. triangle), we have three coordinates and one set of data
  // data will be written in triplicate, but no need to store all that
  // connectivity follows naturally from construction
  Vector<Real>         eltArea(nElts);
  Vector<Real>         eltVol(nElts);
  Vector<Vector<dim3>> surfLocs(nElts);
  Vector<Vector<Real>> surfAvg(nElts);
  Vector<Vector<Real>> surfInt(nElts);
  Vector<Vector<Real>> surfDer(nElts);
  int iElt,iCorner,d;
  Vector<int3> sIdx(nElts);
  // the surface is the mid point of the stream
  int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero    
  // set locations

  int dumpPKZstreams = 0; pp.query("dumpPKZstreams", dumpPKZstreams);
  Vector<Vector<Vector<Real>>> PKZstreams;
  Vector<Real> areaZoneStreams;
  Vector<Real> lsStreams;
  Vector<Real> sVect;
  if(dumpPKZstreams) {
    PKZstreams.resize(nComps-AMREX_SPACEDIM+1);
    lsStreams.resize(6);
    for (int i=0; i<nComps-AMREX_SPACEDIM+1; i++) {
      PKZstreams[i].resize(6);
      for (int j=0; j<6; j++) {
        PKZstreams[i][j].resize(nPtsOnStream);
        for (int k=0; k<nPtsOnStream; k++) {
          PKZstreams[i][j][k] = 0.0;
        }
      }
    }
    sVect.resize(nPtsOnStream);
    areaZoneStreams.resize(6);
    lsStreams.resize(6);
    for (int i=0; i<6;i++) {
      areaZoneStreams[i] = 0;
      lsStreams[i] = 0.0;
    }
  }

  
  Print() << "Setting locations, making triangles, resizing arrays ..." << std::endl;
  
  // get the three stream indices from connectivity faceData
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(iElt,iCorner,d) shared(surfLocs,sIdx,surfAvg,surfInt,surfDer)
#endif    
  for (iElt=0; iElt<nElts; iElt++) {
    surfLocs[iElt].resize(AMREX_SPACEDIM);
    surfAvg[iElt].resize(nAvg);
    surfInt[iElt].resize(nInt);
    surfDer[iElt].resize(nDerOut);
    // define the spatial location of the three corners of the triangle
    for (iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) { 
      sIdx[iElt][iCorner] = faceData[iElt*AMREX_SPACEDIM+iCorner];
      for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	surfLocs[iElt][iCorner][d] = streamData[sIdx[iElt][iCorner]][nPtsOnStream*d+surfPt];
      }
    }
  }

  // evaluate area
  Print() << "Evaluating areas and volumes..." << std::endl;
  Real surfaceArea=0.;
  dim3 A,B,C;
  
  dim3 D,E,F;
  int iPt;
  Real totalVol=0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(iElt,iCorner,A,B,C,D,E,F,d,iPt) shared(eltArea,eltVol,streamData) reduction(+: surfaceArea,totalVol)
#endif
  for (iElt=0; iElt<nElts; iElt++) {
    // set up triangle ABC
    for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
      A[d] = streamData[sIdx[iElt][0]][nPtsOnStream*d+surfPt];
      B[d] = streamData[sIdx[iElt][1]][nPtsOnStream*d+surfPt];
      C[d] = streamData[sIdx[iElt][2]][nPtsOnStream*d+surfPt];
    }
    // find area
    eltArea[iElt] = wedge_area(A,B,C);
    // keep a running total
    surfaceArea+=eltArea[iElt];
    
    eltVol[iElt]=0.;
    for (iPt=1; iPt<nPtsOnStream; iPt++) {
      for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	A[d] = streamData[sIdx[iElt][0]][nPtsOnStream*d+iPt-1];
	B[d] = streamData[sIdx[iElt][1]][nPtsOnStream*d+iPt-1];
	C[d] = streamData[sIdx[iElt][2]][nPtsOnStream*d+iPt-1];
	D[d] = streamData[sIdx[iElt][0]][nPtsOnStream*d+iPt];
	E[d] = streamData[sIdx[iElt][1]][nPtsOnStream*d+iPt];
	F[d] = streamData[sIdx[iElt][2]][nPtsOnStream*d+iPt];
      }
      eltVol[iElt] += wedge_volume(A,B,C,D,E,F);
      totalVol += wedge_volume(A,B,C,D,E,F);
    }
  }
  Print() << "   ... total surface area = " << surfaceArea << std::endl;
  Print() << " ... total volume = " << totalVol << std::endl;
    
  

  // calculate averages
  //Print() << "Calculating averages ..." << std::endl;
  int iAvg, iInt, iDerFlag;
  int3 localSIdx;
  Vector<Vector<Real>> localStreamData;
  Real areaLoc;
  Real filels = 0;
  Real filess = 0;
  Real fileEbar = 0;
  int getEbar = 1; pp.query("getEbar", getEbar);
  int strainIdx = -1;
  if(getEbar) {
    strainIdx = getVarIdx("StrainRate_prog_"+fuelName,avgComps,nAvg);
    amrex::Print() << "Getting Ebar, index: " << strainIdx << std::endl;
  }
  
  
  Print() << "Iterating over elements ..." << std::endl;
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+: filess,filels,fileEbar) private(iElt,d,localSIdx,localStreamData,iAvg,iInt,iDerFlag,areaLoc) shared(surfAvg,nAvg,surfInt,nInt,surfDer,nDerFlag,streamData,nPtsOnStream,derIdxIn,derIdxOut,derCompsIn,pkzFF)
#endif
  for (iElt=0; iElt<nElts; iElt++) {
    //get thread-local stream index, data and element area
    localSIdx = sIdx[iElt];
    localStreamData.resize(AMREX_SPACEDIM);
    for (d=0; d<AMREX_SPACEDIM; d++) {
      localStreamData[d] = streamData[localSIdx[d]];
    }
    areaLoc = eltArea[iElt];

    // evaluate average of each component over the three points
    for (iAvg=0; iAvg<nAvg; iAvg++) {
      surfAvg[iElt][iAvg] = calcAvgVal(avgIdx[iAvg],nPtsOnStream,localStreamData);
    }
    if (getEbar) {
      fileEbar += surfAvg[iElt][strainIdx]*areaLoc;
    }
    for (iInt=0; iInt<nInt; iInt++) {
      surfInt[iElt][iInt] = calcIntegral(intIdx[iInt],nPtsOnStream,localStreamData,areaLoc);
    }
    int outComp;
    for (iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
     if (derCompsIn[iDerFlag]=="flameThickness") {
       outComp = derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp] = deltaT/calcMax(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData);
	// sum for mean
       filels += surfDer[iElt][outComp] * areaLoc;
     }
     if (derCompsIn[iDerFlag]=="flameSpeed") {
       outComp = derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp] = calcIntegral(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData,areaLoc)/rhoY;
       // sum for mean
       filess += surfDer[iElt][outComp] * areaLoc;
     }
     if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
       Real avgMk = calcAvgVal(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData);
       Real avgGk = calcAvgVal(derIdxIn[iDerFlag][1],nPtsOnStream,localStreamData);
       Real det = sqrt(fabs(avgMk*avgMk-avgGk));
       Real k1 = avgMk + det;
       Real k2 = avgMk - det;
       Real zone = 0;
       if   ( sqrt(k1*k1+k2*k2) <  pkzFF         )  zone = 1; // FF
       else {
	 if (                k2 >  half*k1       )  zone = 2; // LP
	 if (       fabs(k2/k1) <= half          )  zone = 3; // LE
	 if (        (-2*k1<k2) && (k2<-half*k1) )  zone = 4; // SP
	 if (       fabs(k1/k2) <= half          )  zone = 5; // TE
	 if (                k1 <  half*k2       )  zone = 6; // TP
       }
       if (zone == 0) {
	 Print() << "Could not find zone" << std::endl;
       }
       outComp = derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp]=k1;
       outComp = derIdxOut[iDerFlag][1];
       surfDer[iElt][outComp]=k2;
       outComp = derIdxOut[iDerFlag][2];
       surfDer[iElt][outComp]=zone;
     }
     
    }
    
  }
  filels /= surfaceArea;
  filess /= surfaceArea;
  fileEbar /= surfaceArea;
  //array reducing stuff (JPDFs, conditionally averaged streamlines)
  //stuff for zone averaging

  int zoneComp, thermalThicknessComp;
  for (iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
    if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
      zoneComp = derIdxOut[iDerFlag][2];
    }
    if (derCompsIn[iDerFlag]=="thermalThickness") {
      thermalThicknessComp = derIdxOut[iDerFlag][0];
    }
  }
  Real ds;
#ifdef _OPENMP
  #pragma omp parallel
  {
#endif //declare local arrays here
    Vector<Vector<Vector<Real>>> PKZstreams_local(PKZstreams);
    Vector<Real> areaZoneStreams_local(areaZoneStreams);
    Vector<Real> lsStreams_local(lsStreams);
#ifdef _OPENMP
#pragma omp for schedule(static) private(iElt,d,localSIdx,localStreamData)
#endif    
    for (iElt=0; iElt<nElts; iElt++) { 
      localSIdx = sIdx[iElt];
      localStreamData.resize(AMREX_SPACEDIM);
      for (d=0; d<AMREX_SPACEDIM; d++) {
	localStreamData[d] = streamData[localSIdx[d]];
      }
      
      areaLoc = eltArea[iElt];
      if (dumpPKZstreams) {
	int zone = surfDer[iElt][zoneComp];
	areaZoneStreams_local[zone-1] += areaLoc; 
	Real thermalThickness = surfDer[iElt][thermalThicknessComp];
	lsStreams[zone-1] += thermalThicknessComp*areaLoc;
	Real surfPoint[3];
	for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++) {
	  Real surfFaceVal = 0.0;
	  for (d=0;d<AMREX_SPACEDIM;d++) {
	    surfFaceVal += localStreamData[d][nPtsOnStream*iComp + surfPt];
	  }
	  surfFaceVal /= (Real)AMREX_SPACEDIM;
	  surfPoint[iComp] = surfFaceVal;
	}
	for (int iPt = 0; iPt<nPtsOnStream; iPt++) {	  
	  Real facePoint[3];	    
	  for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++) {
	    Real faceVal = 0.0;
	    for (d=0;d<AMREX_SPACEDIM;d++) {
	      faceVal += localStreamData[d][nPtsOnStream*iComp + iPt];
	    }
	    faceVal /= (Real)AMREX_SPACEDIM;
	    facePoint[iComp] = faceVal;	    
	  }
	  Real dist = 0.0;
	  for (d=0;d<AMREX_SPACEDIM;d++) {
	    dist += (facePoint[d]-surfPoint[d])*(facePoint[d]-surfPoint[d]);
	  }
	  dist = std::sqrt(dist);
	  if (iPt > surfPt) {
	    dist *= -1;
	  }
	  PKZstreams_local[0][zone-1][iPt] += dist*areaLoc; 
	}
	//probably swap these loop round but cba
	for (int iComp=AMREX_SPACEDIM; iComp<nComps; iComp++) {
	  for (int iPt = 0; iPt<nPtsOnStream; iPt++) {
	    Real faceVal = 0.0;
	    for (d=0;d<AMREX_SPACEDIM;d++) {
	      faceVal += localStreamData[d][nPtsOnStream*iComp + iPt];
	    }
	    faceVal /= (Real)AMREX_SPACEDIM;
	    PKZstreams_local[iComp-AMREX_SPACEDIM+1][zone-1][iPt] += faceVal*areaLoc;
	  }
	}
      }
    }
#ifdef _OPENMP  
#pragma omp critical
  {
#endif
    if (dumpPKZstreams) {
      
      for (int z=0;z<6;z++) {
	areaZoneStreams[z] += areaZoneStreams_local[z];
	lsStreams[z] += lsStreams_local[z];
	for (int iComp=0; iComp<nComps-AMREX_SPACEDIM+1; iComp++) {
	  for (int iPt = 0; iPt<nPtsOnStream; iPt++) {
	    PKZstreams[iComp][z][iPt] += PKZstreams_local[iComp][z][iPt];
	  }
	}
      }
    }
#ifdef _OPENMP          
  }
  }
#endif
  
  //dump characteristic values for this file
  std::string filename=infile+"/characteristics.dat";

  std::ofstream os(filename.c_str(),std::ios::out);
  
  os << filels << " " << filess << " " << fileEbar << std::endl;

  os.close();


  if (dumpPKZstreams) {
    std::string filename;
    Vector<std::string> zoneNames = {"FF","LP","LE","SP","TE","TP"};        

    for (int z = 0; z<6; z++) {
      std::ofstream lsos(filename.c_str(),std::ios::out);
      filename = infile+"/"+zoneNames[z]+"_ls.dat";
      if (areaZoneStreams[z] > 0) {
	lsStreams[z] /= areaZoneStreams[z];
      }
      lsos << lsStreams[z];
      lsos.close();

      filename = infile+"/"+zoneNames[z]+"_distfromseed.dat";
      std::ofstream dsos(filename.c_str(),std::ios::out);
      for (int iPt = nPtsOnStream-1; iPt>=0; iPt--) {
	if (areaZoneStreams[z] > 0) {
	  PKZstreams[0][z][iPt] /= areaZoneStreams[z];
	}
	dsos << PKZstreams[0][z][iPt] << " ";
      }
      dsos.close();
	
      for (int n = AMREX_SPACEDIM; n<nComps; n++) {
	filename = infile+"/"+zoneNames[z]+"_"+variableNames[n]+".dat";
	std::ofstream os(filename.c_str(),std::ios::out);
        for (int iPt = nPtsOnStream-1; iPt>=0; iPt--) {
	  if (areaZoneStreams[z] > 0) {
	    PKZstreams[n-AMREX_SPACEDIM+1][z][iPt] /= areaZoneStreams[z];
	  }
          os << PKZstreams[n-AMREX_SPACEDIM+1][z][iPt] << " ";
        }
        os.close();
      }
    }
  }

  
  // write surface
  Print() << "Writing surface ..." << std::endl;
  writeSurfaceTecplot(infile,
		      nElts, eltArea,  eltVol,  surfLocs,
		      nAvg,  avgComps, surfAvg,
		      nInt,  intComps, surfInt,
		      nDerOut,  derCompsOut, surfDer);
  Print() << "   ... done" << std::endl;

  // write basic
  if (writeBasic) {
    Print() << "Writing basic file ..." << std::endl;
    writeSurfaceBasic(infile,
		      nElts, eltArea,  eltVol,  surfLocs,
		      nAvg,  avgComps, surfAvg,
		      nInt,  intComps, surfInt,
		      nDerOut,  derCompsOut, surfDer);
    Print() << "   ... done" << std::endl;
  }

  return(0);
}

int getVarIdx(std::string varName, Vector<std::string> variableNames, int nComps) {  
    int varIdx=-1;
    for (int iComp=0; iComp<nComps; iComp++)
      if (variableNames[iComp]==varName) varIdx=iComp;
    if (varIdx==-1) {
      std::string msg=varName+" not in stream file";
      Abort(msg);
    }
    return varIdx;
}

Real calcAvgVal(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData)  {
      Real avgVal = 0.0;
      int surfPt = (nPtsOnStream-1)/2;
      // evaluate average of each component over the three points
      for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	avgVal += streamData[iCorner][nPtsOnStream*compIdx+surfPt];
      }
      avgVal *= third;
      return avgVal;
}

Real calcMax(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData) {
  
      Real maxVal = 0.0;
      for (int iPt=1; iPt<nPtsOnStream; iPt++) {
	Real avgVal = 0.0;
	for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	  avgVal += streamData[iCorner][nPtsOnStream*compIdx+iPt];
	}
	avgVal *= third;
	maxVal = max(maxVal,avgVal);
      }
      return maxVal;
}

Real calcIntegral(int compIdx, int nPtsOnStream, Vector<Vector<Real>> streamData, Real eltArea) {

  Real integral = 0.;
  // integrate
  dim3 A,B,C,D,E,F;
  for (int iPt=1; iPt<nPtsOnStream; iPt++) {
    for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
      A[d] = streamData[0][nPtsOnStream*d+iPt-1];
      B[d] = streamData[1][nPtsOnStream*d+iPt-1];
      C[d] = streamData[2][nPtsOnStream*d+iPt-1];
      D[d] = streamData[0][nPtsOnStream*d+iPt];
      E[d] = streamData[1][nPtsOnStream*d+iPt];
      F[d] = streamData[2][nPtsOnStream*d+iPt];
    }
    
    Real vA = streamData[0][nPtsOnStream*compIdx+iPt-1];
    Real vB = streamData[1][nPtsOnStream*compIdx+iPt-1];
    Real vC = streamData[2][nPtsOnStream*compIdx+iPt-1];
    Real vD = streamData[0][nPtsOnStream*compIdx+iPt];
    Real vE = streamData[1][nPtsOnStream*compIdx+iPt];
    Real vF = streamData[2][nPtsOnStream*compIdx+iPt];
    integral += wedge_volume_int(A,vA,B,vB,C,vC,D,vD,E,vE,F,vF);
  } 
  integral /= eltArea;
  return integral;
}



//
// break up the variable list
//
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return amrex::Tokenize(line,std::string(" "));
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

//
// write the stream data to matlab
//
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

//
// write the stream data (midpoint surface) straight to a tecplot file
//
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

//
// write all the surface quantities to a tecplot file
//
void
writeSurfaceTecplot(std::string infile,
		    int& nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		    Vector<Vector<dim3>>& surfLocs,
		    int& nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		    int& nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		    int& nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  std::string vars("VARIABLES = X Y Z area volume");
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
  os << std::setprecision(12);
  for (int iElt=0; iElt<nElts; iElt++) {
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
      // coordinate
      for (int d=0; d<AMREX_SPACEDIM; d++)
	os << surfLocs[iElt][iCorner][d] << " ";
      // area
      os << eltArea[iElt] << " ";
      // volume
      os << eltVol[iElt] << " ";
      // averages
      for (int iAvg=0; iAvg<nAvg; iAvg++)
	os << surfAvg[iElt][iAvg] << " ";
      // integrals
      for (int iInt=0; iInt<nInt; iInt++)
	os << surfInt[iElt][iInt] << " ";
      // derived
      for (int iDer=0; iDer<nDer; iDer++)
	os << surfDer[iElt][iDer] << " ";
      os << std::endl;
    }
  }

  // write connectivity
  for (int iElt=1; iElt<3*nElts;)
    os << iElt++ << " " << iElt++ << " " << iElt++ << std::endl;

  os.close();
}

//
// write all the surface quantities to a tecplot file
//
void
writeSurfaceBasic(std::string infile,
		  int& nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		  Vector<Vector<dim3>>& surfLocs,
		  int& nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		  int& nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		  int& nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  // write averages
  os << std::setprecision(12);
  for (int iElt=0; iElt<nElts; iElt++) {
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
      // coordinate
      for (int d=0; d<AMREX_SPACEDIM; d++)
	os << surfLocs[iElt][iCorner][d] << " ";
      // area
      os << eltArea[iElt] << " ";
      // volume
      os << eltVol[iElt] << " ";
      // averages
      for (int iAvg=0; iAvg<nAvg; iAvg++)
	os << surfAvg[iElt][iAvg] << " ";
      // integrals
      for (int iInt=0; iInt<nInt; iInt++)
	os << surfInt[iElt][iInt] << " ";
      // derived
      for (int iDer=0; iDer<nDer; iDer++)
	os << surfDer[iElt][iDer] << " ";
      os << std::endl;
    }
  }

  os.close();
}

//
// area of the triangle
//
Real
wedge_area(const dim3& A,
	   const dim3& B,
	   const dim3& C)
{
  Real result;  
  Real R1[3], R2[3], R3[3];
  
  for (int i=0; i<3; ++i) {
    R1[i] = B[i] - A[i];
    R2[i] = C[i] - A[i];
  }
  
  R3[0] = R1[1]*R2[2] - R2[1]*R1[2];
  R3[1] = R1[2]*R2[0] - R2[2]*R1[0];
  R3[2] = R1[0]*R2[1] - R2[0]*R1[1];

  result=0;
  for (int i=0; i<3; ++i) {
    result += R3[i]*R3[i];
  }
  result = half*std::sqrt(result);

  return result;
}

//
// volume of a tetrahedron
//
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

    return std::abs(res*sixth);
}

//
// volume of the irregular triangular prism
//
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

//
// integral over the irregular triangular prism
//
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
		       (vB+vC+vD+vF)*vol_BCDF )*quarter;
  
  const Real int_2 = ( (vD+vA+vB+vC)*vol_DABC + 
		       (vC+vD+vE+vF)*vol_CDEF + 
		       (vB+vC+vD+vE)*vol_BCDE )*quarter;
  
  const Real int_3 = ( (vE+vA+vB+vC)*vol_EABC + 
		       (vA+vD+vE+vF)*vol_ADEF + 
		       (vA+vC+vE+vF)*vol_ACEF )*quarter;
  
  const Real int_4 = ( (vE+vA+vB+vC)*vol_EABC + 
		       (vC+vD+vE+vF)*vol_CDEF + 
		       (vA+vC+vE+vD)*vol_ACED )*quarter;
  
  const Real int_5 = ( (vF+vA+vB+vC)*vol_FABC + 
		       (vA+vD+vE+vF)*vol_ADEF + 
		       (vA+vB+vE+vF)*vol_ABEF )*quarter;
  
  const Real int_6 = ( (vF+vA+vB+vC)*vol_FABC + 
		       (vB+vD+vE+vF)*vol_BDEF + 
		       (vA+vB+vD+vF)*vol_ABDF )*quarter;
  
  result = (int_1 + int_2 + int_3 + int_4 + int_5 + int_6)*sixth;
  
  return result;
}

