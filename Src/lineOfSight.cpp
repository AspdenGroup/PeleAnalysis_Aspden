#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

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
  //Vector<std::string> plotVarNames = amrData.PlotVarNames();
  int nVars= pp.countval("vars");
  Vector<std::string> vars(nVars);
  pp.getarr("vars", vars);
  Vector<int> destFillComps(nVars);
  for (int n = 0; n<nVars; n++) {
    destFillComps[n] = n;
  }

  int dir = 0;
  pp.query("dir",dir);
  int dir1 = (dir+1)%3;
  int dir2 = (dir+2)%3;
    
  //MultiFab outdata;
  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel", finestLevel);
  int Nlev = finestLevel + 1;
  Vector<MultiFab*> indata(Nlev);
  
  //structure for outdata
  Box probDomain = amrData.ProbDomain()[finestLevel];
  Real length = amrData.ProbSize()[dir];
  int ldir1 = probDomain.length(dir1);
  int ldir2 = probDomain.length(dir2);
  Vector<Real> tmp(ldir2,0);
  Vector<Vector<Real>> tmp2(ldir1,tmp);
  Vector<Vector<Vector<Real>>> outdata(nVars,tmp2);
  for (int lev = 0; lev < Nlev; lev++) {
    BoxArray probBoxArray = amrData.boxArray(lev);
    indata[lev] = new MultiFab(probBoxArray,DistributionMapping(probBoxArray),nVars+1,0);
    indata[lev]->setVal(1.0,nVars,1);
  }
  Print() << "indata MF allocated" << std::endl;
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
  for (int lev = 0; lev < Nlev; lev++) {
    Real dx = amrData.DxLevel()[lev][0];
    Print() << "Loading data on level " << lev << std::endl;
    amrData.FillVar(*indata[lev],lev,vars,destFillComps);
    Print() << "Data loaded" << std::endl;
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev]); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
	  if (inbox(i,j,k,nVars) > 1e-8) {
	    Vector<int> d(3);
	    d[0] = i;
	    d[1] = j;
	    d[2] = k;
	    for (int n = 0; n < nVars; n++) {
	      outdata[n][d[dir1]][d[dir2]] += dx*inbox(i,j,k,n);
	    }
	  }
	});
    }
  }
  int IO = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealSum(outdata[0][0].data(),nVars*ldir1*ldir2,IO);    
    
  Print() << "Integration completed" << std::endl;
  Vector<Real> d1(ldir1);
  Vector<Real> d2(ldir2);
  Vector<Real> plo = amrData.ProbLo();
  Vector<Real> phi = amrData.ProbHi();
  Real dx = amrData.DxLevel()[finestLevel][0];
  for (int i = 0; i < ldir1; i++) {
    d1[i] = plo[dir1] + (i+0.5)*dx;
  }
  for (int i = 0; i < ldir2; i++) {
    d2[i] = plo[dir2] + (i+0.5)*dx;
  }
  
  if (ParallelDescriptor::IOProcessor())
    {
      std::string outfile=infile+"_LOS_d1_dir"+std::to_string(dir); 
      FILE *file = fopen(outfile.c_str(),"w");
      for (int i = 0; i < ldir1; i++) {
	fprintf(file,"%e ",d1[i]);
      }
      fclose(file);
      outfile=infile+"_LOS_d2_dir"+std::to_string(dir);
      file = fopen(outfile.c_str(),"w");
      for (int i = 0; i < ldir2; i++) {
	fprintf(file,"%e ",d2[i]);
      }
      fclose(file);
      for (int n = 0; n < nVars; n++) {
	outfile=infile+"_LOS_"+vars[n]+"_dir"+std::to_string(dir);
	file = fopen(outfile.c_str(),"w");
	Print() << "outdata[" << n << "].size() = " << outdata[n].size() << std::endl; 
	for (int i = 0; i < ldir1; i++) {
	  Print() << "outdata[" << n << "][" << i << "].size = " << outdata[n][i].size() << std::endl;
	  for (int j = 0; j < ldir2; j++) {
	    fprintf(file,"%e ",outdata[n][i][j]/length);
	  }
	  fprintf(file, "\n");
	}
	fclose(file); 
      }
    }
  
}
















