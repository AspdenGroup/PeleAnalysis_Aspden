#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

void integrate1d(int dir, int dir1, int dir2, int ldir1, int ldir2, Vector<Vector<Vector<Real>>>& outdata, Vector<Real> x, Vector<Real> y, AmrData& amrData, Vector<MultiFab*> indata, int nVars, int finestLevel, int cComp, Real cMin, Real cMax) {
  Real length = 0.0;
  Vector<int> d(3);
  for (int lev = 0; lev <= finestLevel; lev++) {
    Real dxLev = amrData.DxLevel()[lev][dir];
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
	  if (inbox(i,j,k,nVars) > 1e-8 && (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
	    length += dxLev;
	    d[0] = i;
	    d[1] = j;
	    d[2] = k;
	    for (int n = 0; n < nVars; n++) {
	      outdata[n][d[dir1]][d[dir2]] += dxLev*inbox(i,j,k,n);
	    }
	  }
	});
    }
  }
  ParallelDescriptor::ReduceRealSum(outdata[0][0].data(),nVars*ldir1*ldir2);
  ParallelDescriptor::ReduceRealSum(length);
  
  for (int n = 0; n<nVars; n++) {
    for (int i = 0; i < ldir1; i++) {
      for (int j = 0; j < ldir2; j++) {
	outdata[n][i][j] /= length;
      }
    }
  }
  Vector<Real> plo = amrData.ProbLo();
  Vector<Real> phi = amrData.ProbHi();
  Real dxFine = amrData.DxLevel()[finestLevel][dir1];
  Real dyFine = amrData.DxLevel()[finestLevel][dir2];
  for (int i = 0; i < ldir1; i++) {
    x[i] = plo[dir1] + (i+0.5)*dxFine;
  }
  for (int i = 0; i < ldir2; i++) {
    y[i] = plo[dir2] + (i+0.5)*dyFine;
  }  
  return;
}

void integrate2d(int dir, int dir1, int dir2, int ldir, Vector<Vector<Real>>& outdata, Vector<Real> x, AmrData& amrData, Vector<MultiFab*> indata,int nVars, int finestLevel, int cComp, Real cMin, Real cMax) {
  Real area = 0.0;
  Vector<int> d(3);
  for (int lev = 0; lev <= finestLevel; lev++) {
    Real dxLev = amrData.DxLevel()[lev][dir1];
    Real dyLev = amrData.DxLevel()[lev][dir2];
    Real areaLev = dxLev*dyLev;
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
	  if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
	    area += areaLev;
	    d[0] = i;
	    d[1] = j;
	    d[2] = k;
	    for (int n = 0; n < nVars; n++) {
	      outdata[n][d[dir]] += areaLev*inbox(i,j,k,n);
	    }
	  }
	});
    }
  }
  ParallelDescriptor::ReduceRealSum(outdata[0].data(),nVars*ldir);
  ParallelDescriptor::ReduceRealSum(area);
  for (int n = 0; n<nVars; n++) {
    for (int i = 0; i < ldir; i++) {
      outdata[n][i] /= area;
    }
  }
  Vector<Real> plo = amrData.ProbLo();
  Vector<Real> phi = amrData.ProbHi();
  Real dxFine = amrData.DxLevel()[finestLevel][dir];
  for (int i = 0; i < ldir; i++) {
    x[i] = plo[dir] + (i+0.5)*dxFine;
  }
  return;
}
  
void integrate3d(Vector<Real>& outdata, AmrData& amrData, Vector<MultiFab*> indata, int nVars, int finestLevel, int cComp, Real cMin, Real cMax) {
  Real volume = 0.0;
  for (int lev = 0; lev <= finestLevel; lev++) {
    Real dxLev = amrData.DxLevel()[lev][0];
    Real dyLev = amrData.DxLevel()[lev][1];
    Real dzLev = amrData.DxLevel()[lev][2];
    Real volLev = dxLev*dyLev*dzLev;
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
	  if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
	    volume += volLev;
	    for (int n = 0; n < nVars; n++) {
	      outdata[n] += volLev*inbox(i,j,k,n);
	    }
	  }
	});
    }
  }
  ParallelDescriptor::ReduceRealSum(outdata.data(),nVars);
  ParallelDescriptor::ReduceRealSum(volume);
  for (int n = 0; n<nVars; n++) {
    outdata[n] /= volume;
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
  unsigned char *buff=(unsigned char*)malloc(3*dim2*dim2*sizeof(char));
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      int bc   = ((dim1-i-1)*dim2+j)*3;
      Real val = vect[i][j];
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
int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);

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
  for (int n = 0; n<nVars; n++) {
    destFillComps[n] = n;
  }
  
  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel", finestLevel);
  int Nlev = finestLevel + 1;
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
    BoxArray baf = (*indata[lev]).boxArray();
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
  Print() << "Intersect determined" << std::endl;
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
  }
  int integralDimension;
  pp.get("integralDimension",integralDimension);
  int dir, dir1, dir2;
  std::string format;
  std::string outfile= infile+"_integral";
  //pp.query("outfile",outfile);
  switch(integralDimension) {
  case 1:
    {
      pp.get("dir",dir);
      dir1 = (dir+1)%3;
      dir2 = (dir+2)%3;
      pp.get("format",format);
      AMREX_ALWAYS_ASSERT(format=="ppm" || format=="dat");
      Box probDomain = amrData.ProbDomain()[finestLevel];
      int ldir1 = probDomain.length(dir1);
      int ldir2 = probDomain.length(dir2);
      Vector<Real> x(ldir1);
      Vector<Real> y(ldir2);
      Vector<Real> tmp1(ldir2,0.0);
      Vector<Vector<Real>> tmp2(ldir1,tmp1);
      Vector<Vector<Vector<Real>>> outdata(nVars,tmp2);
      //do 1d integration
      integrate1d(dir,dir1,dir2,ldir1,ldir2,outdata,x,y,amrData,indata,nVars,finestLevel,cComp,cMin,cMax);
      Print() << "Integration completed" << std::endl;
      //output data in desired format
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(x,outfile+"_x.dat",ldir1);
	  writeDat1D(y,outfile+"_y.dat",ldir2);
	  for (int n = 0; n < nVars; n++) {
	    writeDat2D(outdata[n],outfile+"_"+vars[n]+".dat",ldir1,ldir2);
	  }
	} else if (format == "ppm") {
	  int goPastMax = 1;
	  pp.query("goPastMax",goPastMax);
	  Vector<Real> vMin(nVars);
	  Vector<Real> vMax(nVars);
	  for (int n=0; n<nVars; n++) {
            char argName[12];
            sprintf(argName,"useminmax%i",n+1);
            int nMinMax = pp.countval(argName);
            if (nMinMax != 2) {
              Abort("Need to specify 2 values for useMinMax");
            } else {
              pp.get(argName, vMin[n], 0);
              pp.get(argName, vMax[n], 1);
            }
	  }
	  for (int n = 0; n < nVars; n++) {
	    writePPM(outdata[n],outfile+"_"+vars[n]+".ppm",ldir1,ldir2,goPastMax,vMin[n],vMax[n]);
	  }
	} else {
	  Abort("Format not recognised!");
	}
      }
      break;
    }
  case 2:
    {
      format="dat"; //probably add an option for binary output
      pp.get("dir1",dir1);
      pp.get("dir2",dir2);
      dir = 3-dir1-dir2;
      Box probDomain = amrData.ProbDomain()[finestLevel];
      int ldir = probDomain.length(dir);
      Vector<Real> x(ldir);
      Vector<Real> tmp(ldir,0.0);
      Vector<Vector<Real>> outdata(nVars,tmp);
      integrate2d(dir,dir1,dir2,ldir,outdata,x,amrData,indata,nVars,finestLevel,cComp,cMin,cMax);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(x,outfile+"_x.dat",ldir);
	  for (int n = 0; n < nVars; n++) {
	    writeDat1D(outdata[n],outfile+"_"+vars[n]+".dat",ldir);
	  }
	} else {
	  Abort("Format not recognised!");
	}
      }
      break;
    }
  case 3:
    {
      format="dat"; //probably add an option for binary output
      Vector<Real> outdata(nVars,0.0);
      integrate3d(outdata,amrData,indata,nVars,finestLevel,cComp,cMin,cMax);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(outdata,outfile+"_allVars.dat",nVars);
	} else {
	  Abort("Format not recognised!");
	}
      }
      break;
    }
  default:
    Abort("integral dimension too high!");
  }
  
  Finalize();
  return 0;  
}
















