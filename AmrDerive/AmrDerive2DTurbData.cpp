#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <new>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <Box.H>
#include <DataServices.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <VisMF.H>

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);
  
  int nProcs = ParallelDescriptor::NProcs();
  if (nProcs>1)
    BoxLib::Error("Not implemented in parallel!");
  int myProc = ParallelDescriptor::MyProc();
  int IOProc = ParallelDescriptor::IOProcessorNumber();
  
  ParmParse pp;
  
  int verbose=0;
  pp.query("verbose",verbose);
  
  int nPlotFiles(pp.countval("infile"));
  if(nPlotFiles <= 0)
    {
      std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
      std::cerr << "Exiting." << std::endl;
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
  
  Array<std::string> plotFileNames(nPlotFiles);
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot)
    pp.get("infile", plotFileNames[iPlot], iPlot);
  
  int nVars(2);
  Array<std::string> varName(nVars);
  varName[0]="x_velocity";
  varName[1]="y_velocity";
  
  Array<int> destFills(nVars);
  for (int v=0; v<nVars; v++ )
    destFills[v] = v;
  
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  FILE *file = fopen("integral_length.dat","w");
  
  double l = 0;
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
    {
      
      std::string infile = plotFileNames[iPlot];
      
      DataServices *dataServices = new DataServices(infile, fileType);
      if( ! dataServices->AmrDataOk())
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
      
      AmrData amrData(dataServices->AmrDataRef());
      
      Real Time      = amrData.Time();
      int  timeSteps = amrData.LevelSteps()[0];
      
      int finestLevel = amrData.FinestLevel();
      
      if (finestLevel!=0)
	{
	  std::cout << "WARNING: Overriding to finest level zero" << std::endl;
	  finestLevel=0;
	}
      
      Box probDomain = amrData.ProbDomain()[finestLevel];
      int ix = probDomain.length(Amrvis::XDIR);
      int jx = probDomain.length(Amrvis::YDIR);
      
      int nCells = ix*jx;

           
      BoxArray domainBoxArray(1);
      domainBoxArray.set(0,probDomain);
      MultiFab mf;
      mf.define(domainBoxArray, nVars, 0, Fab_allocate);
      
      Array<Real> probLo=amrData.ProbLo();
      Array<Real> probHi=amrData.ProbHi();
      
      Real Lx   = probHi[0]-probLo[0];
      Real Ly   = probHi[1]-probLo[1];
      Real dx   = Lx/(Real)ix;
      Real dy   = Ly/(Real)jx;  
      Real dxyz = dx*dy;
      
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Loading plotfile..." << std::endl;
      
      amrData.FillVar(mf, finestLevel, varName, destFills);
      
      for (int n=0; n<nVars; n++)
	amrData.FlushGrids(amrData.StateNumber(varName[n]));
      
      
      // Writes the var as a 1D array
      double *ux_data_array = mf[0].dataPtr(0);
      double *uy_data_array = mf[0].dataPtr(1);
      
      // Calculate u_bar
      double uxb = 0.;
      double uyb = 0.;
      double ux2 = 0.;
      double uy2 = 0.;
      for (int cell=0; cell<nCells; cell++) 
	{
	  double ux = ux_data_array[cell];
	  double uy = uy_data_array[cell];
	  uxb += ux;
	  uyb += uy;
	  ux2 += ux*ux;
	  uy2 += uy*uy;
	}
      // average of velocity
      uxb /= (double)nCells;
      uyb /= (double)nCells;
      // average of the square of velocity
      ux2 /= (double)nCells;
      uy2 /= (double)nCells;
      // rms velocity fluctuation
      double u_bar_x = uxb;
      double u_bar_y = uyb;
      double ux_rms = sqrt( ux2 - uxb*uxb );
      double uy_rms = sqrt( uy2 - uyb*uyb );

      double k = 0.5 * (ux2 + uy2);

      

      
      
      std::cout << "u bar x = "<< u_bar_x << std::endl;
      std::cout << "u bar y = "<< u_bar_y << std::endl;
      std::cout << "ux rms = "<< ux_rms << std::endl;
      std::cout << "uy rms = "<< uy_rms << std::endl;
      std::cout << "k = "<< k << std::endl;
      
      // Looping over entire domain for all r values
      std::cout << "Looping for f(x)" << std::endl;
      double *Qxx = (double*)malloc(ix*sizeof(double));
      double *rx  = (double*)malloc(ix*sizeof(double));
      double Qlxx, Qmxx;
      for (int r=0; r<ix; r++)
	{
	  Qmxx = 0;
	  for (int i=0; i<ix; i++)
	    {
	      Qlxx = 0;
	      int ipr = (i + r) % ix;
	      for (int j=0; j<jx; j++)
		{
		  int cell = (j*ix) + i;
		  int celldx = (j*ix) + ipr;
		  
		  double uxa = ux_data_array[cell];
		  double uxb = ux_data_array[celldx] - u_bar_x;   
		  
		  Qlxx += uxa * uxb;
		}
	      Qmxx += Qlxx / ((double) jx);
	    }
	  Qxx[r] = Qmxx / ((double) ix);             // Mean Qxx[r]
	  Qxx[r] = Qxx[r] / pow(ux_rms,2);           // Calculates fx[r]
	  rx[r] = r*dx;    	  
	}
      int plot_number = iPlot;
      std::string plot_number_str=std::to_string(plot_number);
      FILE *fx_file;
      std::string filename="fx_"+plot_number_str+".dat";
      
      fx_file = fopen(filename.c_str(),"w");
      for (int r=0; r<ix; r++)
	{
	  fprintf(fx_file,"%e %e \n",Qxx[r],rx[r]);
	}
      fclose(fx_file);
      
      std::cout << "Looping for f(y)" << std::endl;
      double *Qyy = (double*)malloc(jx*sizeof(double));
      double *ry  = (double*)malloc(jx*sizeof(double));
      double Qlyy, Qmyy;
      
      for (int r=0; r<jx; r++)
	{
	  Qmyy = 0;
	  for (int i=0; i<ix; i++)
	    {
	      Qlyy = 0;
	      for (int j=0; j<jx; j++)
		{
		  int jpr = (j + r) % jx;
     
		  int cell = (j*ix) + i;
		  int celldy = (jpr*ix) + i;
		  
		  double uya = uy_data_array[cell];
		  double uyb = uy_data_array[celldy] - u_bar_y;
		  
		  Qlyy += uya * uyb;
		}
	      Qmyy += Qlyy / ((double) jx);
	    }
	  Qyy[r] = Qmyy / ((double) ix);           // Mean Qyy[r]
	  Qyy[r] = Qyy[r] / pow(uy_rms,2);         // Calculates fy[r]
	  ry[r] = r*dy;
	}
      
      // Trapezium Rule
      std::cout << "Running Trapezium Rule" << std::endl;
      double sum;
      sum = 0.;
      for (int r=1; r<ix-1; r++)
	{
	  sum += Qxx[r];
	}
      double lx = 0.5 * dx * ((Qxx[0] + Qxx[ix-1]) + 2 * sum);
      sum = 0.;
      for (int r=1; r<jx-1; r++)
	{
	  sum += Qyy[r];
	}
      double ly = 0.5 * dy * ((Qyy[0] + Qyy[ix-1]) + 2 * sum);

      fprintf(file,"%e %e %e %e %e %e \n",Time,lx,ly,(lx+ly)/2,ux_rms, k);     
      std::cout << "integral length x = "<< lx << std::endl;
      std::cout << "integral length y = "<< ly << std::endl;
    }
  fclose(file);
  BoxLib::Finalize();    
  return(0);
}










