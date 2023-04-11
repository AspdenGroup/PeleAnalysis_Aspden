#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

void minTheta(int dirx, int diry, int dirz, Vector<Vector<Vector<Real>>>& outdata, Vector<Real>& r, Vector<Real>& z, AmrData& amrData, Vector<MultiFab*> indata, int nVars, int finestLevel, int cComp, Real cMin, Real cMax) {

  Box probDomain = amrData.ProbDomain()[finestLevel];
  int ldirr = (int)(std::sqrt(2)*std::max(probDomain.length(dirx),probDomain.length(diry))/2.0);
  int ldirz = probDomain.length(dirz);
  int bandingFactor = 1;
  int refRatio = 1;
  Vector<Real> plo = amrData.ProbLo();
  for (int lev = finestLevel; lev >= 0; lev--) {
    Real dzLev = amrData.DxLevel()[lev][dirz];
    Real dxLev = amrData.DxLevel()[lev][dirx];
    Real dyLev = amrData.DxLevel()[lev][diry];
    if (lev < finestLevel) refRatio *= amrData.RefRatio()[lev];
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
	  if (inbox(i,j,k,nVars) > 1e-8 && (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {	    
	    Real xlo = plo[0] + i*dxLev;
	    Real ylo = plo[1] + j*dyLev;
	    for (int bx = 0; bx < bandingFactor; bx++) {
	      for(int by = 0; by < bandingFactor; by++) {
		Real x = xlo + ((bx+0.5)*dxLev)/(Real)bandingFactor;
		Real y = ylo + ((by+0.5)*dyLev)/(Real)bandingFactor;
		Real r = std::sqrt(x*x + y*y);	     
		int ridx = (int)(r/dxLev);
		if (ridx < ldirr) {
		  for (int rx = 0; rx < refRatio; rx++) {
		    for (int ry = 0; ry < refRatio; ry++) {
		      for (int n = 0; n < nVars; n++) {
			outdata[n][refRatio*ridx+rx][refRatio*k+ry] = std::min(outdata[n][refRatio*ridx+rx][refRatio*k+ry],inbox(i,j,k,n-1));
		      }
		    }
		  }
		}
	      }
	    }
	  });
	
	}
	
    }
  }
  for (int i = 0; i < ldirr; i++) {
    for (int n = 0; n<nVars+1; n++) {
      ParallelDescriptor::ReduceRealMin(outdata[n][i].data(),ldirz);
    }
  }
  Real dxFine = amrData.DxLevel()[finestLevel][dirx];
  Real dzFine = amrData.DxLevel()[finestLevel][dirz];
  for (int i = 0; i < ldirr; i++) {
    r[i] = (i+0.5)*dxFine;
  }
  for (int i = 0; i < ldirz; i++) {
    z[i] = plo[dirz] + (i+0.5)*dzFine;
  }  
  return;
}

void writeDat1D(Vector<Real> vect, std::string filename, int dim) {
  FILE *file = fopen(filename.c_str(),"w");
  for (int i = 0; i < dim; i++) {
    fprintf(file,"%e ",vect[i]);
  }
  fclose(file);
  return;
}

void writeDat2D(Vector<Vector<Real>> vect, std::string filename, int dim1, int dim2) {
  FILE *file = fopen(filename.c_str(),"w");
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      fprintf(file,"%e ",vect[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return;
}

void writePPM(Vector<Vector<Real>> vect, std::string filename, int dim1, int dim2, int goPastMax, Real vMin, Real vMax) {
  unsigned char *buff=(unsigned char*)malloc(3*dim1*dim2*sizeof(char));
  for (int iz=0; iz<dim2; iz++) {
    for (int ir=0; ir<dim1; ir++) {
      int bc   = ((dim2-iz-1)*dim1+ir)*3;
      Real val = vect[ir][iz];
      Real colour = fmax(0.,fmin(1.5,(val-vMin)/(vMax-vMin)));
      if (colour<0.125) {
	buff[bc]   = 0;
	buff[bc+1] = 0;
	buff[bc+2] = (int)((colour+0.125)*1020.);
      } else if (colour<0.375)  {
	buff[bc]   = 0;
	buff[bc+1] = (int)((colour-0.125)*1020.);
	buff[bc+2] = 255;
      } else if (colour<0.625)  {
	buff[bc]   = (int)((colour-0.375)*1020.);
	buff[bc+1] = 255;
	buff[bc+2] = (int)((0.625-colour)*1020.);
      } else if (colour<0.875)  {
	buff[bc]   = 255;
	buff[bc+1] = (int)((0.875-colour)*1020.);
	buff[bc+2] = 0;
      } else if (colour<1.000)  {
	buff[bc]   = (int)((1.125-colour)*1020.);
	buff[bc+1] = 0;
	buff[bc+2] = 0;
      } else if (goPastMax==1) {
	if (colour<1.125)  {
	  buff[bc]   = (int)((colour-0.875)*1020.);
	  buff[bc+1] = 0;
	  buff[bc+2] = (int)((colour-1.000)*1020.);
	} else if (colour<1.250) {
	  buff[bc]   = 255;
	  buff[bc+1] = 0;
	  buff[bc+2] = (int)((colour-1.000)*1020.);
	} else if (colour<1.500)  {
	  buff[bc]   = 255;
	  buff[bc+1] = (int)((colour-1.250)*1020.);
	  buff[bc+2] = 255;
	} else { // default if above 1.5 with goPastMax==1
	  buff[bc]   = 255;
	  buff[bc+1] = 255;
	  buff[bc+2] = 255;
	}
      } else { // default if above 1 with goPastMax==0
	buff[bc]   = 128;
	buff[bc+1] = 0;
	buff[bc+2] = 0;
      }
    }
  }
  FILE *file = fopen(filename.c_str(),"w");
  fprintf(file,"P6\n%i %i\n255\n",dim1,dim2);
  fwrite(buff,dim1*dim2*3,sizeof(unsigned char),file);
  fclose(file);
  return;
}

void findMinMax(Vector<Vector<Real>> vect, int dim1, int dim2, Real &min, Real &max) {
  min = vect[0][0];
  max = vect[0][0];
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      if (vect[i][j] < min) min = vect[i][j];
      if (vect[i][j] > max) max = vect[i][j];
    }
  }
  return;
}

int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);
  {
    ParmParse pp;
  
    std::string infile;
    pp.get("infile",infile);
    Print() << "infile = " << infile << std::endl; 
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();
    
    // read in the variable names to use
    int nVars= pp.countval("vars");
    Vector<std::string> vars(nVars);
    pp.getarr("vars", vars);
    Vector<int> destFillComps(nVars);
    Print() << "nVars= " << nVars << std::endl;
    for (int n = 0; n<nVars; n++) {
      destFillComps[n] = n;
      Print() << "var[" << n << "]= " << vars[n] << std::endl;
    }
  
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;
    std::string cVar;
    Real cMin,cMax;
    int cComp=-1;
    pp.query("cVar",cVar);
    pp.query("cMin",cMin);
    pp.query("cMax",cMax);
    if (!cVar.empty()) {
      for (int n = 0; n<nVars; n++) {
	if (vars[n] == cVar) {
	  cComp = n;
	  break;
	}
      }
      if (cComp < 0) {
	Abort("cVar not in list of vars!");
      }
    }
    int dir, dir1, dir2;
    std::string format="dat";
    pp.query("format",format);
    
    std::string outfile= infile+"_min_rz";
    if(!cVar.empty()) {
      outfile+="_c"+cVar+"_"+std::to_string(cMin)+"_"+std::to_string(cMax);
    }    

    Vector<MultiFab*> indata(Nlev);
    for (int lev = 0; lev < Nlev; lev++) {
      BoxArray probBoxArray = amrData.boxArray(lev);
      indata[lev] = new MultiFab(probBoxArray,DistributionMapping(probBoxArray),nVars+1,0);
      Print() << "Loading data on level " << lev << std::endl;
      amrData.FillVar(*indata[lev],lev,vars,destFillComps);
      Print() << "Data loaded" << std::endl;
      indata[lev]->setVal(1.0,nVars,1);
    }
    Print() << "Determining intersects..." << std::endl;
    for (int lev = 0; lev < finestLevel; lev++) {
      BoxArray baf = (*indata[lev+1]).boxArray();
      baf.coarsen(amrData.RefRatio()[lev]);	  
      for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
	FArrayBox& myFab = (*indata[lev])[mfi];
	int idx = mfi.index();      
	std::vector< std::pair<int,Box> > isects = baf.intersections((*indata[lev]).boxArray()[idx]);
	for (int ii = 0; ii < isects.size(); ii++) {
	  myFab.setVal(0.0,isects[ii].second,nVars,1);
	}
      }
    }
    Print() << "Intersects determined" << std::endl;
    int dirx = 0;
    int diry = 1;
    int dirz = 2;
    
    Box probDomain = amrData.ProbDomain()[finestLevel];
    int ldirr = (int)(std::sqrt(2)*std::max(probDomain.length(dirx),probDomain.length(diry))/2.0);
    Print() << "Length of domain in r: " << ldirr << std::endl;
    int ldirz = probDomain.length(dirz);
    Vector<Real> r(ldirr);
    Vector<Real> z(ldirz);
    Vector<Real> tmp1(ldirz,0.0);
    Vector<Vector<Real>> tmp2(ldirr,tmp1);
    Vector<Vector<Vector<Real>>> outdata(nVars,tmp2);
    //do min
    minTheta(dirx,diry,dirz,outdata,r,z,amrData,indata,nVars,finestLevel,cComp,cMin,cMax);
    Print() << "Minimum check completed" << std::endl;
    //output data in desired format
    Print() << "Writing data as "+format << std::endl;
    if (ParallelDescriptor::IOProcessor()) {
      if (format == "dat") {
	writeDat1D(r,outfile+"_r.dat",ldirr);
	writeDat1D(z,outfile+"_z.dat",ldirz);
	for (int n = 0; n < nVars; n++) {
	  writeDat2D(outdata[n],outfile+"_"+vars[n-1]+".dat",ldirr,ldirz);
	}
      } else if (format == "ppm") {
	int goPastMax = 1;
	pp.query("goPastMax",goPastMax);
	Vector<Real> vMin(nVars);
	Vector<Real> vMax(nVars);
	findMinMax(outdata[0],ldirr,ldirz,vMin[0],vMax[0]);
	for (int n=0; n<nVars; n++) {
	  char argName[12];
	  sprintf(argName,"useminmax%i",n+1);
	  int nMinMax = pp.countval(argName);
	  if (nMinMax > 0) {
	    Print() << "Reading min/max from command line" << std::endl;
	    if (nMinMax != 2) {
	      Abort("Need to specify 2 values for useMinMax");
	    } else {
	      pp.get(argName, vMin[n], 0);
	      pp.get(argName, vMax[n], 1);
	    }
	  } else {
	    Print() << "Using file values for min/max" << std::endl;
	    findMinMax(outdata[n],ldirr,ldirz,vMin[n],vMax[n]);
	  }
	}
	for (int n = 0; n < nVars; n++) { 
	  writePPM(outdata[n],outfile+"_"+vars[n]+".ppm",ldirr,ldirz,goPastMax,vMin[n],vMax[n]);
	}
      } //can add more formats here if we want - add to assert above
    }
  } 
  Finalize();
  return 0;  
}

















