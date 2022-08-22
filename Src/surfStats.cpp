#include <string>
#include <iostream>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

#include <surfStats.H>

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  if (ParallelDescriptor::NProcs()>1)
    Abort("Code is not yet parallel safe");
  
  ParmParse pp;

  // read infile name from inputs
  Vector<std::string> infiles;
  int nFiles = pp.countval("infiles");
  infiles.resize(nFiles);
  pp.getarr("infiles",infiles);
  Print () << "nFiles= " << nFiles << std::endl;
  //make arrays to hold data
  Real ls = 0.0;
  Real ss = 0.0;
  Real Ebar = 0.0;
  int nBin = 100;
  Real minL,maxL,minS,maxS,minK,maxK,minE,maxE,mink1,maxk1,mink2,maxk2,minStretch,maxStretch;
  int turb = 0;
  pp.query("turb",turb);
  pp.query("nBin",nBin);
  Real rhoY;
  Real ls_input = -1;
  Real ss_input = -1;
  Real EbarT = -1;
  std::string fuelName = "H2";
  pp.query("fuelName", fuelName);
  std::string strainComp = "StrainRate_prog_"+fuelName;
  pp.query("strainComp",strainComp);
  std::string curvatureComp="MeanCurvature_prog_"+fuelName;
  pp.query("curvatureComp",curvatureComp);
  std::string FCRComp=fuelName+"_ConsumptionRate";
  pp.query("FCRComp",FCRComp);    

  pp.get("rhoY",rhoY);  
  pp.query("ls",ls_input);
  pp.query("ss",ss_input);
  pp.query("Ebar",EbarT);
  if (ls_input < 0 || ss_input < 0 || EbarT < 0) {
    Real omega2 = 0.0;
    int pressure_regime = 0;
    Real sL = 0.0;
    Real lL = 0.0;
    Real Ka = 0.0;
    pp.get("omega2",omega2);
    pp.get("pr",pressure_regime);
    pp.get("sl",sL);
    pp.get("ll",lL);
    if (pressure_regime == 0) {
      ss_input = sL*std::exp(0.08*omega2);
      ls_input = lL*std::exp(-0.06*omega2);
      EbarT = 1.7+0.05*omega2;
    } else {
      ss_input = sL*(1+0.47*omega2);
      ls_input = lL/(1+0.26*omega2);
      EbarT = 1+0.13*omega2;
    }    
    if (turb > 0) {
      pp.get("ka",Ka);
      ss_input *= 1+0.26*std::exp(-0.038*omega2)*std::sqrt(Ka);
      ls_input *= 1+0.22*std::exp(-0.026*omega2)*std::sqrt(Ka);
    }
    
  }
		
  
  Real M,M1,M2;
  pp.get("M",M);
  pp.get("M1",M1);
  pp.get("M2",M2);
  
  
  Real tau_input = ls_input/ss_input;
  Real Ebar_input = EbarT;// /tau_input;
  Real reacTemp, prodTemp;
  pp.get("reacTemp",reacTemp);
  pp.get("prodTemp",prodTemp);
  Real const deltaT=prodTemp-reacTemp;
  std::string pkzMkComp = curvatureComp;
  std::string pkzGkComp = "GaussianCurvature_prog_"+fuelName;
  pp.query("pkzMkVar",pkzMkComp);
  pp.query("pkzGkVar",pkzGkComp);    
  std::string gradTComp = "ModGradTemp";
  pp.query("gradTComp",gradTComp);
  
  if (turb == 0) {
    minL = 0.0;
    maxL = 1.5;
    minS = 0.0;
    maxS = 5;
    minK = -1.5;
    maxK = 1.5;
    minE = -3;
    maxE = 12;
    mink1 = minK;
    maxk1 = maxK;
    mink2 = minK;
    maxk2 = maxK;
    minStretch = -1;
    maxStretch = 3;
  } else {
    minL = 0.0;
    maxL = 1.5;
    minS = 0.0;
    maxS = 5;
    minK = -1.5;
    maxK = 1.5;
    minE = -3;
    maxE = 12;
    mink1 = minK;
    maxk1 = maxK;
    mink2 = minK;
    maxk2 = maxK;
    minStretch = -1;
    maxStretch = 3;
  }
  
    
  pp.query("minL",minL);
  minL *= ls_input;
  pp.query("maxL",maxL);
  maxL *= ls_input;
  pp.query("minS",minS);
  minS *= ss_input;
  pp.query("maxS",maxS);
  maxS *= ss_input;
  pp.query("minK",minK);
  minK /= ls_input;
  pp.query("maxK",maxK);
  maxK /= ls_input; 
  pp.query("minE",minE);
  minE /= tau_input;
  pp.query("maxE",maxE);
  maxE /= tau_input;
  pp.query("mink1",mink1);
  mink1 /= ls_input;
  pp.query("maxk1",maxk1);
  maxk1 /= ls_input;
  pp.query("mink2",mink2);
  mink2 /= ls_input;
  pp.query("maxk2",maxk2);
  maxk2 /= ls_input;
  pp.query("maxStretch",maxStretch);
  pp.query("minStretch",minStretch);
  
  Vector<Real> sBins(nBin);
  Real ds = (maxS-minS)/(Real)nBin; 
  Vector<Real> lBins(nBin);
  Real dl = (maxL-minL)/(Real)nBin; 
  Vector<Real> kBins(nBin);
  Real dk = (maxK-minK)/(Real)nBin; 
  Vector<Real> EBins(nBin);
  Real de = (maxE-minE)/(Real)nBin; 
  Vector<Real> k1Bins(nBin);
  Real dk1 = (maxk1-mink1)/(Real)nBin;
  Vector<Real> k2Bins(nBin);
  Real dk2 = (maxk2-mink2)/(Real)nBin; 
  Vector<Real> stretchBins(nBin);
  Real dS = (maxStretch-minStretch)/(Real)nBin;
  
  Vector<Real> fractionalContribution(12);
  for (int iBin = 0; iBin < nBin; iBin++) {
    sBins[iBin] = minS + (Real)iBin*ds;    
    lBins[iBin] = minL + (Real)iBin*dl;
    kBins[iBin] = minK + (Real)iBin*dk;
    EBins[iBin] = minE + (Real)iBin*de;
    k1Bins[iBin] = mink1 + (Real)iBin*dk1;
    k2Bins[iBin] = mink2 + (Real)iBin*dk2;
    stretchBins[iBin] = minStretch + (Real)iBin*dS;
  }
  Vector<Vector<Real>> ksJPDF(nBin);
  Vector<Vector<Real>> EsJPDF(nBin); 
  Vector<Vector<Real>> k1k2JPDF(nBin);
  Vector<Vector<Real>> k1k2FMJPDF(nBin);
  Vector<Vector<Vector<Real>>> kEsJPDF(nBin);
  Vector<Vector<Real>> stretchsJPDF(nBin);
  Vector<Vector<Real>> stretchIndsJPDF(nBin);
  
  for (int iBin = 0; iBin < nBin; iBin++) {
    
    ksJPDF[iBin].resize(nBin);
    EsJPDF[iBin].resize(nBin);
    k1k2JPDF[iBin].resize(nBin);
    k1k2FMJPDF[iBin].resize(nBin);
    kEsJPDF[iBin].resize(nBin);
    stretchsJPDF[iBin].resize(nBin);
    stretchIndsJPDF[iBin].resize(nBin);
    for (int j = 0; j < nBin; j++) {
      ksJPDF[iBin][j] = 0.0;
      EsJPDF[iBin][j] = 0.0;
      k1k2JPDF[iBin][j] = 0.0;
      k1k2FMJPDF[iBin][j] = 0.0;
      stretchsJPDF[iBin][j] = 0.0;
      stretchIndsJPDF[iBin][j] = 0.0;
      kEsJPDF[iBin][j].resize(nBin);
      for (int k = 0; k<nBin; k++) {
	kEsJPDF[iBin][j][k] = 0.0;
      }
    }
  }
  for (int i=0;i<12;i++) {
    fractionalContribution[i] = 0.0;
  }
  

  Real pkzLength = ls_input;
  pp.query("pkzLength",pkzLength);
  Real pkzFF = half/pkzLength;
  
  Real filels, filess,fileEbar;
  for (int iFile=0; iFile<nFiles; iFile++) {
    filels = 0.0;
    filess = 0.0;
    fileEbar = 0.0;
    Print() << "Processing " << infiles[iFile] << "..." << std::endl; 
    //declare stuff
    int nStreams, nElts, nPtsOnStream, nComps;
    Vector<std::string> variableNames;
    Vector<int> faceData;
    Vector<Vector<Real>> streamData;
  
    // Read file
    readStreamBin(infiles[iFile], nStreams, nElts, nPtsOnStream, nComps,
		  variableNames, faceData, streamData);
    // report
    Print() << "nStreams      = " << nStreams << std::endl;
    Print() << "nElements     = " << nElts << std::endl;
    Print() << "nPtsOnStreams = " << nPtsOnStream << std::endl;
    int strainIdx = getVarIdx(strainComp,variableNames,nComps);
    int curvatureIdx = getVarIdx(curvatureComp,variableNames,nComps);
    int FCRIdx = getVarIdx(FCRComp,variableNames,nComps);
    int gradTIdx = getVarIdx(gradTComp,variableNames,nComps);
    int pkzMkIdx = getVarIdx(pkzMkComp,variableNames,nComps);
    int pkzGkIdx = getVarIdx(pkzGkComp,variableNames,nComps);
    // make space to hold output surface
    // for each element (i.e. triangle), we have three coordinates and one set of data
    // data will be written in triplicate, but no need to store all that
    // connectivity follows naturally from construction
    Vector<Real>         eltArea(nElts);
    Vector<Real>         eltVol(nElts);
    Vector<Vector<dim3>> surfLocs(nElts);

    int iElt,iCorner;
    Vector<int3> sIdx(nElts);
    int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero    
    // set locations
    Print() << "Setting locations, making triangles etc ..." << std::endl;
    int idxStream, d;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(iElt,idxStream,iCorner,d) shared(surfLocs,sIdx)
#endif    
    for (iElt=0; iElt<nElts; iElt++) {
      surfLocs[iElt].resize(AMREX_SPACEDIM);
      for (iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	sIdx[iElt][iCorner] = faceData[iElt*AMREX_SPACEDIM+iCorner];
	idxStream=faceData[iElt*AMREX_SPACEDIM+iCorner];
	for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	  surfLocs[iElt][iCorner][d] = streamData[idxStream][nPtsOnStream*d+surfPt];
	}
      }
    }
  

    // evaluate area
    Print() << "Evaluating areas and volumes..." << std::endl;
    Real surfaceArea=0.;
    int3 iStream;
    dim3 A,B,C;

    dim3 D,E,F;
    int iPt;
    Real totalVol=0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(iElt,iCorner,iStream,A,B,C,D,E,F,d,iPt) shared(eltArea,eltVol) reduction(+: surfaceArea,totalVol)
#endif
    for (iElt=0; iElt<nElts; iElt++) {
      for (iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
	iStream[iCorner]=sIdx[iElt][iCorner];
      for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	A[d] = streamData[iStream[0]][nPtsOnStream*d+surfPt];
	B[d] = streamData[iStream[1]][nPtsOnStream*d+surfPt];
	C[d] = streamData[iStream[2]][nPtsOnStream*d+surfPt];
      }
      // find area
      eltArea[iElt] = wedge_area(A,B,C);
      // keep a running total
      surfaceArea+=eltArea[iElt];

      eltVol[iElt]=0.;
      for (iPt=1; iPt<nPtsOnStream; iPt++) {
	for (d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	  A[d] = streamData[iStream[0]][nPtsOnStream*d+iPt-1];
	  B[d] = streamData[iStream[1]][nPtsOnStream*d+iPt-1];
	  C[d] = streamData[iStream[2]][nPtsOnStream*d+iPt-1];
	  D[d] = streamData[iStream[0]][nPtsOnStream*d+iPt];
	  E[d] = streamData[iStream[1]][nPtsOnStream*d+iPt];
	  F[d] = streamData[iStream[2]][nPtsOnStream*d+iPt];
	}
	eltVol[iElt] += wedge_volume(A,B,C,D,E,F);
	totalVol += wedge_volume(A,B,C,D,E,F);
      }
    }
    Print() << "   ... total surface area = " << surfaceArea << std::endl;
    Print() << " ... total volume = " << totalVol << std::endl;
    
    Print() << "Iterating over elements ..." << std::endl;


#ifdef _OPENMP
    
#pragma omp parallel
{
#endif
  Real curvatureLoc, strainLoc, speedLoc, thermalThicknessLoc,maxModGradT, avgMk,avgGk,det,k1,k2, contr,areaLoc,stretchLoc, stretchIndLoc;
  int curvatureBinIdx, strainBinIdx, speedBinIdx, stretchBinIdx, stretchIndBinIdx ,k1BinIdx, k2BinIdx, zone, zoneIdx;
  int3 localSIdx;
  Vector<Vector<Real>> localStreamData;
  localStreamData.resize(AMREX_SPACEDIM);
  
  Vector<Vector<Real>> ksJPDF_file(nBin);
  Vector<Vector<Real>> EsJPDF_file(nBin); 
  Vector<Vector<Real>> k1k2JPDF_file(nBin);
  Vector<Vector<Real>> k1k2FMJPDF_file(nBin);
  Vector<Vector<Vector<Real>>> kEsJPDF_file(nBin);
  Vector<Vector<Real>> stretchsJPDF_file(nBin);
  Vector<Vector<Real>> stretchIndsJPDF_file(nBin);
  Vector<Real> fractionalContribution_file(12);

  
  for (int iBin = 0; iBin < nBin; iBin++) {
    
    ksJPDF_file[iBin].resize(nBin);
    EsJPDF_file[iBin].resize(nBin);
    k1k2JPDF_file[iBin].resize(nBin);
    k1k2FMJPDF_file[iBin].resize(nBin);
    kEsJPDF_file[iBin].resize(nBin);
    stretchsJPDF_file[iBin].resize(nBin);
    stretchIndsJPDF_file[iBin].resize(nBin);
    for (int j = 0; j < nBin; j++) {
      ksJPDF_file[iBin][j] = 0.0;
      EsJPDF_file[iBin][j] = 0.0;
      k1k2JPDF_file[iBin][j] = 0.0;
      k1k2FMJPDF_file[iBin][j] = 0.0;
      stretchsJPDF_file[iBin][j] = 0.0;
      stretchIndsJPDF_file[iBin][j] = 0.0;
      kEsJPDF_file[iBin][j].resize(nBin);
      for (int k = 0; k<nBin; k++) {
	kEsJPDF_file[iBin][j][k] = 0.0;
      }
    }
  }
  for (int i=0;i<12;i++) {
    fractionalContribution_file[i] = 0.0;
  }
  
#ifdef _OPENMP
#pragma omp for schedule(static) reduction(+: filess,filels) private(iElt,d)
#endif
    for (iElt=0; iElt<nElts; iElt++) {
      localSIdx = sIdx[iElt];
      for (d=0; d<AMREX_SPACEDIM; d++) {
	localStreamData[d] = streamData[localSIdx[d]];
      }
      areaLoc = eltArea[iElt];
      strainLoc = calcAvgVal(strainIdx,nPtsOnStream,localStreamData);
      strainBinIdx = floor((strainLoc-minE)/de);
      fileEbar += strainLoc*areaLoc/surfaceArea;
      
      curvatureLoc = calcAvgVal(curvatureIdx,nPtsOnStream,localStreamData);
      curvatureBinIdx = floor((curvatureLoc-minK)/dk);

      stretchLoc = M*(curvatureLoc*ls_input + (strainLoc-Ebar_input)*tau_input);
      stretchBinIdx = floor((stretchLoc-minStretch)/dS);
      stretchIndLoc = M1*curvatureLoc*ls_input + M2*(strainLoc-Ebar_input)*tau_input;
      stretchIndBinIdx = floor((stretchIndLoc-minStretch)/dS);
      speedLoc = calcIntegral(FCRIdx,nPtsOnStream,localStreamData,areaLoc);
      speedLoc /= rhoY;
      if (std::isnan(speedLoc)) {
	continue;
      }
      speedBinIdx = floor((speedLoc-minS)/ds);
      filess += speedLoc*areaLoc/surfaceArea;
      
      
      maxModGradT = calcMax(gradTIdx,nPtsOnStream,localStreamData);
      thermalThicknessLoc = deltaT/maxModGradT;
      filels += thermalThicknessLoc * areaLoc/surfaceArea;
      
      //Print() << filels << std::endl;
      // also need to know which variable to use as mean and gaussian curvature
      avgMk = calcAvgVal(pkzMkIdx,nPtsOnStream,localStreamData);
      avgGk = calcAvgVal(pkzGkIdx,nPtsOnStream,localStreamData);
      
      det = sqrt(fabs(avgMk*avgMk-avgGk));
      k1 = avgMk + det;
      k2 = avgMk - det;
      k1BinIdx = (int)((k1-mink1)/dk1);
      k2BinIdx = (int)((k2-mink2)/dk2); 
      zone = 0;
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
      zoneIdx = zone-1;
      
      contr = areaLoc/surfaceArea;
      if (speedBinIdx >= 0 && speedBinIdx < nBin) {
	if (curvatureBinIdx >=  0 && curvatureBinIdx < nBin) {
	  ksJPDF_file[curvatureBinIdx][speedBinIdx] += contr; 
	  if (strainBinIdx >= 0 && strainBinIdx < nBin) {
	    kEsJPDF_file[curvatureBinIdx][strainBinIdx][speedBinIdx] += contr;
	  }
	}
	if (strainBinIdx >= 0 && strainBinIdx < nBin) {
	  EsJPDF_file[strainBinIdx][speedBinIdx] += contr;
	}
	if (stretchBinIdx >= 0 && stretchBinIdx < nBin) {
	  stretchsJPDF_file[stretchBinIdx][speedBinIdx] += contr; 
	}
	if (stretchIndBinIdx >= 0 && stretchIndBinIdx < nBin) {
	  stretchIndsJPDF_file[stretchIndBinIdx][speedBinIdx] += contr;
	}
      }
      if (k1BinIdx >= 0 && k1BinIdx < nBin && k2BinIdx >=0 && k2BinIdx < nBin) {
	k1k2JPDF_file[k1BinIdx][k2BinIdx] += contr;
	k1k2FMJPDF_file[k1BinIdx][k2BinIdx] += speedLoc*contr;
      }
      fractionalContribution_file[zoneIdx] += contr;      
      fractionalContribution_file[6+zoneIdx] += speedLoc*contr;
    }
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      for (int iBin1 = 0; iBin1 < nBin; iBin1++) {
	for (int iBin2 = 0; iBin2 < nBin; iBin2++) {
	  ksJPDF[iBin1][iBin2] += ksJPDF_file[iBin1][iBin2]/(Real)nFiles;
	  EsJPDF[iBin1][iBin2] += EsJPDF_file[iBin1][iBin2]/(Real)nFiles;
	  k1k2JPDF[iBin1][iBin2] += k1k2JPDF_file[iBin1][iBin2]/(Real)nFiles;
	  k1k2FMJPDF[iBin1][iBin2] += k1k2FMJPDF_file[iBin1][iBin2]/(Real)nFiles;
	  stretchsJPDF[iBin1][iBin2] += stretchsJPDF_file[iBin1][iBin2]/(Real)nFiles;
	  stretchIndsJPDF[iBin1][iBin2] += stretchIndsJPDF_file[iBin1][iBin2]/(Real)nFiles;
	  for (int iBin3 = 0; iBin3 < nBin; iBin3++) {
	    kEsJPDF[iBin1][iBin2][iBin3] += kEsJPDF_file[iBin1][iBin2][iBin3]/(Real)nFiles;
	  }
	}
      }
      for (int i = 0; i<12; i++) {
	if (!std::isnan(fractionalContribution_file[i])) {
	  fractionalContribution[i] += fractionalContribution_file[i]/(Real)nFiles;
	}
      }
#ifdef _OPENMP      
    }
 }
#endif
 
      Print() << "Finished iterating elements" << std::endl;
      ss += filess;
      Print() << "file ss = " << filess << std::endl;
      ls += filels;
      Print() << "file ls = " << filels << std::endl;
      Ebar += fileEbar;
  }
  
  ss /= (Real)nFiles;
  ls /= (Real)nFiles;
  Ebar /= (Real)nFiles;
  Vector<Real> characteristics(3);
  characteristics[0] = ls;
  characteristics[1] = ss;
  characteristics[2] = Ebar*ls/ss;
  checkUnity1(fractionalContribution,6);
  checkUnity2(ksJPDF,nBin,nBin);
  checkUnity2(EsJPDF,nBin,nBin);
  checkUnity2(k1k2JPDF,nBin,nBin);
  checkUnity2(stretchsJPDF,nBin,nBin);
  checkUnity2(stretchIndsJPDF,nBin,nBin);
  //writeSingleArray(characteristics,3,"properties/characteristics");
  writeSingleArray(sBins,nBin,"JPDF/s.dat");
  writeSingleArray(EBins,nBin,"JPDF/e.dat");
  writeSingleArray(kBins,nBin,"JPDF/k.dat");
  //writeSingleArray(stretchBins,nBin,"JPDF/stretch,dat");
  writeSingleArray(fractionalContribution,12,"JPDF/fractionalContribution.dat");
  writeDoubleArray(ksJPDF,nBin,nBin,"JPDF/ksJPDF.dat");
  writeDoubleArray(EsJPDF,nBin,nBin,"JPDF/EsJPDF.dat");
  writeDoubleArray(k1k2JPDF,nBin,nBin,"JPDF/k1k2JPDF.dat");
  writeDoubleArray(k1k2FMJPDF,nBin,nBin,"JPDF/k1k2FMJPDF.dat");
  writeDoubleArray(stretchsJPDF,nBin,nBin,"JPDF/stretchsJPDF.dat");
  writeDoubleArray(stretchIndsJPDF,nBin,nBin,"JPDF/stretchIndsJPDF.dat");
  for (int n = 0; n<nBin; n++) {
    std::string filename = "JPDF/KES/JPDF"+std::to_string(n)+".dat";
    writeDoubleArray(kEsJPDF[n],nBin,nBin,filename);
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
/*
Real* calcArea() {

}


Real* calcVolume() {


}
*/

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
  
  int iFile;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(iFile) shared(streamData) num_threads(nFiles)
#endif
  for (iFile=0; iFile<nFiles; iFile++) {
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

void
checkUnity1(Vector<Real> arr, int arrsize) {
  Real sum = 0;
  for (int i = 0; i < arrsize; i++) {
    sum += arr[i];
  }
  std::cout << "sum = " << sum << std::endl;
  return;
}


void
checkUnity2(Vector<Vector<Real>> arr, int arrsize1, int arrsize2) {
  Real sum = 0;
  for (int i = 0; i < arrsize1; i++) {
    for (int j = 0; j < arrsize2; j++) {
      sum += arr[i][j];
      //Print() << arr[i][j] << std::endl;
    }
  }
  std::cout << "sum = " << sum << std::endl;
  return;
}

void
writeSingleArray(Vector<Real> arr, int arrsize, std::string filename) {
  FILE *file;
  file = fopen(filename.c_str(),"w");
  for (int i = 0; i < arrsize; i++) {
    fprintf(file,"%e ",arr[i]);
  }
  fclose(file);
  return;
}

void
writeDoubleArray(Vector<Vector<Real>> arr, int arrsize1, int arrsize2, std::string filename) {
  FILE *file;
  file = fopen(filename.c_str(),"w");
  for (int j = 0; j < arrsize2; j++) {
    for (int i = 0; i < arrsize1; i++) {
      fprintf(file,"%e ",arr[i][j]);
    }
    fprintf(file,"\n");
  }
  fclose(file);
  return;
}
/*
void
writeTripleArray(Vector<Vector<Vector<Real>>> arr, int arrsize1, int arrsize2, int arrsize3, std::string filename) {
  FILE *file;
  file = fopen(filename.c_str(),"w");
  for (int k = 0; k < arrsize3; k++) {
    for (int j = 0; j < arrsize2; j++) {
      for (int i = 0; i < arrsize1; i++) {
	fprintf(file,"%e ",arr[i][j]);
      }
      fprintf(file,"\n");
    }
  }
  fclose(file);
  return;
}
*/
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

