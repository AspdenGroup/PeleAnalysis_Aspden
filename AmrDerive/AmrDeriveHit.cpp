
//
// This is a version of AmrDerive.cpp that calculates integrals of
// quantities and writes out scalars instead of plotfiles.
//


#include <new>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <unistd.h>

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "VisMF.H"
#include "Geometry.H"

#include "AmrDeriveHit_F.H"

// This one is for the ghost cells (filcc)
#include "xtra_F.H"

#ifdef DO_SPECTRA_PARALLEL
#  define DO_SPECTRA
#  if ( defined(FRANKLIN) || defined(JAGWAR))
#    include "dfftw_mpi.h"
#  else
#    include "fftw_mpi.h"
#  endif
#endif

#ifdef DO_SPECTRA_SERIAL
#define DO_SPECTRA
#include "fftw3.h"
#endif

#ifdef DO_SPECTRA
Real ****mallocCorrltnData(int nPlotFiles, int ix, int jx, int kx);
void freeCorrltnData(Real ****CorrltnData, int nPlotFiles);
#endif

bool verbose=false;

#ifdef DO_SPECTRA

#define MODIFY_SPECTRAundef

static
void
EvalIntLen(Real *l, Real ***Q, int *hx, const Real *dx)
{
  if (verbose)
    std::cout << "   + evaluating integral length scales..." << std::endl;
  for (int c=0; c<DIMS; c++) {                              // loop through each velocity component
    int flag = 1;
    Real fOld = Q[c][c][0]+1.;
    l[c] = 0.;
    for (int i=0; (i<hx[c])&&(flag); i++) {                 // integrate over spatial extent
      Real f = Q[c][c][i];
      if (f<fOld && f>0.) {                                 // but only as far as the first minimum or zero
	l[c] += f*dx[c];
	fOld = f;
      } else {
	flag = 0;
      }
    }
  }
}

static
void
EvalKlmgrv(Real *Klmgrv, Real *spec, int *hx, const Real *dx)
{
    if (verbose)
	std::cout << "   + evaluating Kolmgorov length scales..." << std::endl;    

    int ix, jx, kx;
    ix = hx[0]<<1;    jx = hx[1]<<1;     kx = hx[2]<<1;

    int wavenumbers = evalWavenumbers(ix,jx,kx); // assume cubic for now

    int wnLo = (wavenumbers*3)/8;
    int wnHi = wnLo*2;
    int wnMid= (int)((wnLo+wnHi)/2);
    /*
    printf("Wavenumbers: %i (Lo %i, Hi %i, Mid %i)\n", wavenumbers, wnLo, wnHi, wnMid);
    printf("Energy[30]: %e %e %e\n",spec[90],spec[91],spec[92]);
    */
    for (int c=0; c<DIMS; c++) {

	Real eta = 3*dx[c];
	Real A   = spec[wnMid*3+c] * exp(5.2*((Real)(wnMid))*eta);

	int iter = 0;
	Real res;
	
	do {
	    Real d1  = 0.;
	    Real d2  = 0.;
	    Real Mee = 0.;
	    Real MAA = 0.;
	    Real MeA = 0.;
	    
	    for (int wn=wnLo; wn<wnHi; wn++) {
		Real kappa = (Real) wn;
		Real E     = spec[wn*3+c];
		Real Exp   = exp(-5.2*kappa*eta);
		Real dfde  = -5.2*kappa*A*Exp;
		Real dfdA  = Exp;

		d1  += dfde*(E-A*Exp);
		d2  += dfdA*(E-A*Exp);

		Mee += dfde*dfde;
		MAA += dfdA*dfdA;
		MeA += dfde*dfdA;
	    }

	    Real det  = Mee*MAA - MeA*MeA;
	    
	    Real dEta = (1./det)*( MAA*d1-MeA*d2);
	    Real dA   = (1./det)*(-MeA*d1+Mee*d2);

	    eta += dEta;
	    A   += dA;
		
	    res = (dEta*dEta + dA*dA)/(eta*eta + A*A);
	    iter++;
	    /*
	    if (verbose)
		std::cout << "      + res(" << iter << ") = " << res << std::endl;
	    */

	} while (res>TOL && iter<ITERS);

	if (iter<ITERS) {
	    std::cout << "      + " << iter << " iterations required." << std::endl;
	    Klmgrv[c] = eta;
	} else {
	    std::cout << "      + maximum iterations reached!" << std::endl;
	    Klmgrv[c] =-eta;
	}

    }

    if (verbose) {
	std::cout << "      + eta(i,j,k)    = " << Klmgrv[0]       << ", " << Klmgrv[1]       << ", " << Klmgrv[2]       << std::endl;
	std::cout << "      + eta(i,j,k)/dx = " << Klmgrv[0]/dx[0] << ", " << Klmgrv[1]/dx[1] << ", " << Klmgrv[2]/dx[2] << std::endl;
    }

    for (int c=0; c<DIMS; c++) {
	Real wols = log(spec[wnHi*3+c]/spec[wnLo*3+c]) / (-5.2*(Real)(wnHi-wnLo));
	printf("Without least squares: %e (%e dx)\n",wols,wols/dx[c]);
    }

}

static
void
Spectra(AmrData &amrData, const MultiFab &mf,
#ifdef DO_SPECTRA_PARALLEL
	fftwnd_mpi_plan &plan,  fftwnd_mpi_plan &inverse_plan, 
#else
	fftw_plan &p,  fftw_plan &ip, 
#endif
	fftw_complex *local_xdata, fftw_complex *local_ydata, fftw_complex *local_zdata,
	Real *spec, Real ***Q, Real *u2)
{
#ifdef DO_SPECTRA_SERIAL
  fftw_plan plan, inverse_plan;
#endif
  int finestLevel(amrData.FinestLevel());
  Box probDomain(amrData.ProbDomain()[finestLevel]);
  int dix = probDomain.length(XDIR);
  int djx = probDomain.length(YDIR);
  int dkx = probDomain.length(ZDIR);
  int wavenumbers = evalWavenumbers(dix,djx,dkx);
  int ix, jx, kx;
  const Real *dx = amrData.DxLevel()[finestLevel].dataPtr();

#if 0
  if (verbose) {
    std::cout << "Outside MFIter (START)" << std::endl
	      << "dix: " << dix << std::endl
	      << "djx: " << djx << std::endl
	      << "dkx: " << dkx << std::endl
	      << "wavenumbers: " << wavenumbers << std::endl
	      << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl
	      << "Outside MFIter (END)" << std::endl << std::endl;
  }
#endif

  if (verbose)
    std::cout << "   + calculating energy spectra..." << std::endl;

#ifdef DO_SPECTRA_SERIAL
  if (verbose)
    std::cout << "   + calculating spectra for the first time..." << std::endl;
  for (int k = 0; k < dkx; k++) {
      for (int j = 0; j < djx; j++) {
	  for (int i = 0; i < dix; i++) {
	      int fft_cell = (k*djx + j)*dix + i;
	      local_xdata[fft_cell][0] = 0.;   local_xdata[fft_cell][1] = 0.;
	      local_ydata[fft_cell][0] = 0.;   local_ydata[fft_cell][1] = 0.;
	      local_zdata[fft_cell][0] = 0.;   local_zdata[fft_cell][1] = 0.;
	  }
      }
  }
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_xdata,local_xdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_xdata,local_xdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_ydata,local_ydata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_ydata,local_ydata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_zdata,local_zdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_zdata,local_zdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  if (verbose)
    std::cout << "      done..." << std::endl;
#endif

  for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
    const FArrayBox &myFab = mf[ntmfi];
    
    const int  *dlo = myFab.loVect();
    const int  *dhi = myFab.hiVect();
    const int   dat_ix = dhi[0] - dlo[0] + 1;
    const int   dat_jx = dhi[1] - dlo[1] + 1;
    const int   dat_kx = dhi[2] - dlo[2] + 1;

    const int  *lo = ntmfi.validbox().loVect();
    const int  *hi = ntmfi.validbox().hiVect();
    ;           ix = hi[0] - lo[0] + 1;
    ;           jx = hi[1] - lo[1] + 1;
    ;           kx = hi[2] - lo[2] + 1;

    int ingrow = (dat_ix - ix)>>1;
    int jngrow = (dat_jx - jx)>>1;
    int kngrow = (dat_kx - kx)>>1;

    const Real* ux = myFab.dataPtr(0);
    const Real* uy = myFab.dataPtr(1);
    const Real* uz = myFab.dataPtr(2);

#if 0
    if (ParallelDescriptor::MyProc()==2)
      std::cout << "Inside MFIter (START)" << std::endl
		<< "ix, jx, kx: " << ix << ", " << jx << ", " << kx << std::endl
		<< "dat_ix, dat_jx, dat_kx: " << dat_ix << ", " << dat_jx << ", " << dat_kx << std::endl
		<< "ingrow, jngrow, kngrow: " << ingrow << ", " << jngrow << ", " << kngrow << std::endl
		<< "Inside MFIter (END)" << std::endl;
#endif

    for (int k = 0; k < kx; k++) {
      int dat_k = k+kngrow;
      for (int j = 0; j < jx; j++) {
	int dat_j = j+jngrow;
	for (int i = 0; i < ix; i++) {
	  int dat_i = i+ingrow;
	  int dat_cell = (dat_k*dat_jx + dat_j)*dat_ix + dat_i;
	  int fft_cell = (k*jx + j)*ix + i;
#ifdef DO_SPECTRA_PARALLEL
	  local_xdata[fft_cell].re = ux[dat_cell];   local_xdata[fft_cell].im = 0.;
	  local_ydata[fft_cell].re = uy[dat_cell];   local_ydata[fft_cell].im = 0.;
	  local_zdata[fft_cell].re = uz[dat_cell];   local_zdata[fft_cell].im = 0.;
#else
	  local_xdata[fft_cell][0] = ux[dat_cell];   local_xdata[fft_cell][1] = 0.;
	  local_ydata[fft_cell][0] = uy[dat_cell];   local_ydata[fft_cell][1] = 0.;
	  local_zdata[fft_cell][0] = uz[dat_cell];   local_zdata[fft_cell][1] = 0.;
#endif
	}
      }
    }
    // Let's see how many processors we're working with here
    int nProcs(ParallelDescriptor::NProcs());
    int myProc = ParallelDescriptor::MyProc();

#ifdef DO_SPECTRA_PARALLEL
    fftwnd_mpi(plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);
#else
    plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
    plan = fftw_plan_dft_3d(kx,jx,ix,local_ydata,local_ydata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
    plan = fftw_plan_dft_3d(kx,jx,ix,local_zdata,local_zdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
#endif

    // Integrate for spectra on this grid
    int  wn;
    int  wni, wnj, wnk, wn2;
    Real Rx, Ry, Rz;
    Real Ix, Iy, Iz;
    Real Mx, My, Mz;
    Real div = (Real) (dix*djx*dkx); // Divisor to normalise transform

    for (int k = 0; k < kx; k++) {
      int kp = k + lo[2];
      if (kp<dkx-kp) wnk = kp;
      else           wnk = (dkx-kp);

      for (int j = 0; j < jx; j++) {
	int jp = j + lo[1];
	if (jp<djx-jp) wnj = jp;
	else           wnj = (djx-jp);
	
	for (int i = 0; i < ix; i++) {
	  int ip = i + lo[0];
	  if (ip<dix-ip) wni = ip;
	  else           wni = (dix-ip);
	  
	  wn2 = wni*wni + wnj*wnj + wnk*wnk;
	  wn = (int) (0.5+sqrt((Real)wn2));

	  int cell = (k*jx + j)*ix + i;	    

#if DO_SPECTRA_PARALLEL
	  Rx = local_xdata[cell].re / div; 	    Ix = local_xdata[cell].im / div;	  Mx = (Rx*Rx + Ix*Ix);   // Use div to normalise the transform
	  Ry = local_ydata[cell].re / div; 	    Iy = local_ydata[cell].im / div;	  My = (Ry*Ry + Iy*Iy);
	  Rz = local_zdata[cell].re / div; 	    Iz = local_zdata[cell].im / div;	  Mz = (Rz*Rz + Iz*Iz);
#else
	  Rx = local_xdata[cell][0] / div; 	    Ix = local_xdata[cell][1] / div;	  Mx = (Rx*Rx + Ix*Ix);   // Use div to normalise the transform
	  Ry = local_ydata[cell][0] / div; 	    Iy = local_ydata[cell][1] / div;	  My = (Ry*Ry + Iy*Iy);
	  Rz = local_zdata[cell][0] / div; 	    Iz = local_zdata[cell][1] / div;	  Mz = (Rz*Rz + Iz*Iz);
#endif

	  /* Integrate for energy spectrum */
	  if (wn<wavenumbers) {
	    spec[wn*3  ] += 0.5 * Mx; // Half for energy
	    spec[wn*3+1] += 0.5 * My;
	    spec[wn*3+2] += 0.5 * Mz;
	  }

	  /* Form energy spectrum tensor to invert for correlation tensor */
#ifdef DO_SPECTRA_PARALLEL
	  local_xdata[cell].re = Mx;	  local_xdata[cell].im = 0.;
	  local_ydata[cell].re = My;	  local_ydata[cell].im = 0.;
	  local_zdata[cell].re = Mz;	  local_zdata[cell].im = 0.;
#else
	  local_xdata[cell][0] = Mx;	  local_xdata[cell][1] = 0.;
	  local_ydata[cell][0] = My;	  local_ydata[cell][1] = 0.;
	  local_zdata[cell][0] = Mz;	  local_zdata[cell][1] = 0.;
#endif
	}
      }
    }

    /* Invert spectrum tensors to get correlation tensors */
#ifdef DO_SPECTRA_PARALLEL
    fftwnd_mpi(inverse_plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);
#else
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_ydata,local_ydata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_zdata,local_zdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
#endif

    // Extract correlation functions we're interested in from tensor
    if (lo[2]==0) { // Have to make sure we're on the bottom slab for x and y correlations
      for (int i = 0; i < dix>>1; i++) {
	int cell = i; // j=k=0
#ifdef DO_SPECTRA_PARALLEL
	Q[0][0][i] = local_xdata[cell].re / u2[0];
	Q[1][0][i] = local_ydata[cell].re / u2[1];
	Q[2][0][i] = local_zdata[cell].re / u2[2];
#else
	Q[0][0][i] = local_xdata[cell][0] / u2[0];
	Q[1][0][i] = local_ydata[cell][0] / u2[1];
	Q[2][0][i] = local_zdata[cell][0] / u2[2];
#endif
      }
      for (int j = 0; j < djx>>1; j++) {
	int cell = j*ix; // i=k=0
#ifdef DO_SPECTRA_PARALLEL
	Q[0][1][j] = local_xdata[cell].re / u2[0];
	Q[1][1][j] = local_ydata[cell].re / u2[1];
	Q[2][1][j] = local_zdata[cell].re / u2[2];
#else
	Q[0][1][j] = local_xdata[cell][0] / u2[0];
	Q[1][1][j] = local_ydata[cell][0] / u2[1];
	Q[2][1][j] = local_zdata[cell][0] / u2[2];
#endif
      }
    }
    // Only do the k bits we have on this processor
    for (int k = 0, kp = lo[2]; ( (k<kx) && (kp<(dkx>>1)) ); k++, kp++) {
      int cell = k*jx*ix; // i=j=0
#ifdef DO_SPECTRA_PARALLEL
      Q[0][2][kp] = local_xdata[cell].re / u2[0];
      Q[1][2][kp] = local_ydata[cell].re / u2[1];
      Q[2][2][kp] = local_zdata[cell].re / u2[2];
#else
      Q[0][2][kp] = local_xdata[cell][0] / u2[0];
      Q[1][2][kp] = local_ydata[cell][0] / u2[1];
      Q[2][2][kp] = local_zdata[cell][0] / u2[2];
#endif
    }

  }
  
  // Reduce spectra
  for (int i=0; i<(wavenumbers*DIMS); i++)
    ParallelDescriptor::ReduceRealSum(spec[i]);

  // Reduce correlation functions
  for (int c=0; c<3; c++) {
    for (int i=0; i<(dix>>1); i++) 
      ParallelDescriptor::ReduceRealSum(Q[c][0][i]);
    for (int j=0; j<(djx>>1); j++)
      ParallelDescriptor::ReduceRealSum(Q[c][1][j]);
    for (int k=0; k<(dkx>>1); k++)
      ParallelDescriptor::ReduceRealSum(Q[c][2][k]);
  }

}
#endif


























#ifdef MODIFY_SPECTRA
static
void
ModifySpectra(AmrData &amrData, MultiFab &mf,
#ifdef DO_SPECTRA_PARALLEL
	fftwnd_mpi_plan &plan,  fftwnd_mpi_plan &inverse_plan, 
#else
	fftw_plan &p,  fftw_plan &ip, 
#endif
	fftw_complex *local_xdata, fftw_complex *local_ydata, fftw_complex *local_zdata)
{
#ifdef DO_SPECTRA_SERIAL
  fftw_plan plan, inverse_plan;
#endif
  int finestLevel(amrData.FinestLevel());
  Box probDomain(amrData.ProbDomain()[finestLevel]);
  int dix = probDomain.length(XDIR);
  int djx = probDomain.length(YDIR);
  int dkx = probDomain.length(ZDIR);
  int wavenumbers = evalWavenumbers(dix,djx,dkx);
  int ix, jx, kx;
  const Real *dx = amrData.DxLevel()[finestLevel].dataPtr();

#if 0
  if (verbose) {
    std::cout << "Outside MFIter (START)" << std::endl
	      << "dix: " << dix << std::endl
	      << "djx: " << djx << std::endl
	      << "dkx: " << dkx << std::endl
	      << "wavenumbers: " << wavenumbers << std::endl
	      << "dx: " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl
	      << "Outside MFIter (END)" << std::endl << std::endl;
  }
#endif

  if (verbose)
    std::cout << "   + calculating energy spectra..." << std::endl;

#ifdef DO_SPECTRA_SERIAL
  if (verbose)
    std::cout << "   + calculating spectra for the first time..." << std::endl;
  for (int k = 0; k < dkx; k++) {
      for (int j = 0; j < djx; j++) {
	  for (int i = 0; i < dix; i++) {
	      int fft_cell = (k*djx + j)*dix + i;
	      local_xdata[fft_cell][0] = j;   local_xdata[fft_cell][1] = k;
	      local_ydata[fft_cell][0] = k;   local_ydata[fft_cell][1] = i;
	      local_zdata[fft_cell][0] = i;   local_zdata[fft_cell][1] = j;
	  }
      }
  }
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_xdata,local_xdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_xdata,local_xdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_ydata,local_ydata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_ydata,local_ydata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  plan         = fftw_plan_dft_3d(dkx,djx,dix,local_zdata,local_zdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);            fftw_destroy_plan(plan);
  inverse_plan = fftw_plan_dft_3d(dkx,djx,dix,local_zdata,local_zdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
  if (verbose)
    std::cout << "      done..." << std::endl;
#endif

  for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
    FArrayBox &myFab = mf[ntmfi];
    
    const int  *dlo = myFab.loVect();
    const int  *dhi = myFab.hiVect();
    const int   dat_ix = dhi[0] - dlo[0] + 1;
    const int   dat_jx = dhi[1] - dlo[1] + 1;
    const int   dat_kx = dhi[2] - dlo[2] + 1;

    const int  *lo = ntmfi.validbox().loVect();
    const int  *hi = ntmfi.validbox().hiVect();
    ;           ix = hi[0] - lo[0] + 1;
    ;           jx = hi[1] - lo[1] + 1;
    ;           kx = hi[2] - lo[2] + 1;

    int ingrow = (dat_ix - ix)>>1;
    int jngrow = (dat_jx - jx)>>1;
    int kngrow = (dat_kx - kx)>>1;

    Real* ux = myFab.dataPtr(0);
    Real* uy = myFab.dataPtr(1);
    Real* uz = myFab.dataPtr(2);

    for (int k = 0; k < kx; k++) {
      int dat_k = k+kngrow;
      for (int j = 0; j < jx; j++) {
	int dat_j = j+jngrow;
	for (int i = 0; i < ix; i++) {
	  int dat_i = i+ingrow;
	  int dat_cell = (dat_k*dat_jx + dat_j)*dat_ix + dat_i;
	  int fft_cell = (k*jx + j)*ix + i;
#ifdef DO_SPECTRA_PARALLEL
	  local_xdata[fft_cell].re = ux[dat_cell];   local_xdata[fft_cell].im = 0.;
	  local_ydata[fft_cell].re = uy[dat_cell];   local_ydata[fft_cell].im = 0.;
	  local_zdata[fft_cell].re = uz[dat_cell];   local_zdata[fft_cell].im = 0.;
#else
	  local_xdata[fft_cell][0] = ux[dat_cell];   local_xdata[fft_cell][1] = 0.;
	  local_ydata[fft_cell][0] = uy[dat_cell];   local_ydata[fft_cell][1] = 0.;
	  local_zdata[fft_cell][0] = uz[dat_cell];   local_zdata[fft_cell][1] = 0.;
#endif
	}
      }
    }
#ifdef DO_SPECTRA_PARALLEL
    fftwnd_mpi(plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);
#else
    plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
    plan = fftw_plan_dft_3d(kx,jx,ix,local_ydata,local_ydata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
    plan = fftw_plan_dft_3d(kx,jx,ix,local_zdata,local_zdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
#endif

    // Modify spectrum
    int  wn;
    int  wni, wnj, wnk, wn2;
    //?? Real wni, wnj, wnk, wn2;
    Real div = (Real) (dix*djx*dkx); // Divisor to normalise transform

    std::cout << "Modifying spectrum..." << std::endl;
    std::cout << wavenumbers << std::endl;

    for (int k = 0; k < kx; k++) {
      int kp = k + lo[2];
      if (kp<dkx-kp) wnk = kp;
      else           wnk = (dkx-kp);

      for (int j = 0; j < jx; j++) {
	int jp = j + lo[1];
	if (jp<djx-jp) wnj = jp;
	else           wnj = (djx-jp);
	
	for (int i = 0; i < ix; i++) {
	  int ip = i + lo[0];
	  if (ip<dix-ip) wni = ip;
	  else           wni = (dix-ip);
	  
	  wn2 = wni*wni + wnj*wnj + wnk*wnk;
	  wn = (int) (0.5+sqrt((Real)wn2));

	  int cell = (k*jx + j)*ix + i;	    

	  Real fact = 1.;

	  if (sqrt(wn2)>(0.9*(Real)wavenumbers))
	      fact = 1.+ (sqrt(wn2)-(0.9*(Real)wavenumbers))/(0.1*wavenumbers);
	  
#ifdef DO_SPECTRA_PARALLEL
	  local_xdata[cell].re *= fact;	  local_xdata[cell].im *= fact;
	  local_ydata[cell].re *= fact;	  local_ydata[cell].im *= fact;
	  local_zdata[cell].re *= fact;	  local_zdata[cell].im *= fact;
#else
	  local_xdata[cell][0] *= fact;	  local_xdata[cell][1] *= fact;
	  local_ydata[cell][0] *= fact;	  local_ydata[cell][1] *= fact;
	  local_zdata[cell][0] *= fact;	  local_zdata[cell][1] *= fact;
#endif
	}
      }
    }

    /* Invert spectrum tensors to get correlation tensors */
#ifdef DO_SPECTRA_PARALLEL
    fftwnd_mpi(inverse_plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);
#else
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_ydata,local_ydata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
    inverse_plan = fftw_plan_dft_3d(kx,jx,ix,local_zdata,local_zdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(inverse_plan);    fftw_destroy_plan(inverse_plan);
#endif

    /* Unload */
    for (int k = 0; k < kx; k++) {
      int dat_k = k+kngrow;
      for (int j = 0; j < jx; j++) {
	int dat_j = j+jngrow;
	for (int i = 0; i < ix; i++) {
	  int dat_i = i+ingrow;
	  int dat_cell = (dat_k*dat_jx + dat_j)*dat_ix + dat_i;
	  int fft_cell = (k*jx + j)*ix + i;
#ifdef DO_SPECTRA_PARALLEL
	  ux[dat_cell] = local_xdata[fft_cell].re/div;
	  uy[dat_cell] = local_ydata[fft_cell].re/div;
	  uz[dat_cell] = local_zdata[fft_cell].re/div;
#else
	  ux[dat_cell] = local_xdata[fft_cell][0]/div;
	  uy[dat_cell] = local_ydata[fft_cell][0]/div;
	  uz[dat_cell] = local_zdata[fft_cell][0]/div;
#endif
	}
      }
    }
  }
}
#endif









static
void
Hit(AmrData &amrData, const MultiFab &mf, Real *AverageData, int *PdfData, Real *PdfMin, Real *PdfMax)
{
  int finestLevel(amrData.FinestLevel());
  const Real *dx = amrData.DxLevel()[finestLevel].dataPtr();

  Box probDomain(amrData.ProbDomain()[finestLevel]);
  int dix = probDomain.length(XDIR);
  int djx = probDomain.length(YDIR);
  int dkx = probDomain.length(ZDIR);
  Real cells = (Real)(dix*djx*dkx);

  Array<Real> hitAvgs(BARVARS, 0.);
  Real *avgs = hitAvgs.dataPtr();
  
  if (verbose)
    std::cout << "   + calculating averages..." << std::endl;

  for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
    const FArrayBox &myFab = mf[ntmfi];
    
    const Real *dat = myFab.dataPtr();
    const int  *dlo = myFab.loVect();
    const int  *dhi = myFab.hiVect();
    const int  *lo  = ntmfi.validbox().loVect();
    const int  *hi  = ntmfi.validbox().hiVect();

    FORT_HIT(dat, ARLIM(dlo), ARLIM(dhi), ARLIM(lo), ARLIM(hi), dx, avgs, PdfData, PdfMin, PdfMax);

  }
  
  for (int k=0; k < BARVARS; k++) {
    ParallelDescriptor::ReduceRealSum(avgs[k]);
    AverageData[k] = avgs[k]/cells;
  }
  
  for (int k=0; k < 12*PDFBINS; k++) {
    ParallelDescriptor::ReduceIntSum(PdfData[k]);
  }
  
  if (verbose) {
    std::cout << "ux, uy, uz: " << AverageData[UXBAR] << ", " << AverageData[UYBAR] << ", " << AverageData[UZBAR] << std::endl
	      << "ux2, uy2, uz2: " << AverageData[UX2BAR] << ", " << AverageData[UY2BAR] << ", " << AverageData[UZ2BAR] << std::endl
	      << "Diss: " << AverageData[DBAR] << std::endl
#ifdef DO_FORCING
	      << "Forcing: " << AverageData[FBAR] << std::endl
#endif
	      << std::endl;
  }
}


int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);
    
  if(argc == 1)
    PrintUsage(argv[0]);

  ParmParse pp;

  if(pp.contains("help"))
    PrintUsage(argv[0]);

  FArrayBox::setFormat(FABio::FAB_IEEE_32);

  if(pp.contains("verbose") || (pp.contains("v"))) {
    verbose = true;
    AmrData::SetVerbose(true);
  }

  // Let's set verbose to true, without setting AmrData(verbose)
  if (ParallelDescriptor::IOProcessor())
    verbose = true;
  else
    verbose = false;

  // Let's see how many processors we're working with here
  int nProcs(ParallelDescriptor::NProcs());    VSHOWVAL(verbose, nProcs);
  int myProc = ParallelDescriptor::MyProc();

  // Count number of input plot files
  int nPlotFiles(pp.countval("infile"));
  if(nPlotFiles <= 0) {
    std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
    std::cerr << "Exiting." << std::endl;
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  VSHOWVAL(verbose, nPlotFiles);

  // Make an array of srings containing paths of input plot files
  Array<std::string> plotFileNames(nPlotFiles);
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    pp.get("infile", plotFileNames[iPlot], iPlot);
    VSHOWVAL(verbose, plotFileNames[iPlot]);
  }

  // Get an output file name
  // Hacked out for now
  //  std::string oFile;
  //  pp.get("outfile", oFile);
  //  VSHOWVAL(verbose, oFile);

  // More random initialisation stuff
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  // To do
  //// Handle multiple plot files properly
  //// Check plot domains are the same
  //// Initialise average array before loop
  ////// then pass it to "Inflow"

  Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
  Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
  Array<Real>           time(nPlotFiles);
    
  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Loading " << plotFileNames[iPlot] << std::endl;
      
    dataServicesPtrArray[iPlot] = new DataServices(plotFileNames[iPlot], fileType);               // Populate DataServices array
      
    if( ! dataServicesPtrArray[iPlot]->AmrDataOk())                                               // Check AmrData ok
      DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not
      
    amrDataPtrArray[iPlot] = &(dataServicesPtrArray[iPlot]->AmrDataRef());                        // Populate DataPtrArray

    time[iPlot] = amrDataPtrArray[iPlot]->Time();

    if (verbose) std::cout << "Time = " << time[iPlot] << std::endl;
  }
    
  // Check that all the plot files have the same domain
  // // Check domain is the same size
  // // Check base grid is the same size
  // // Find highest refinement level across plots
  // // Pick a level at which to make the slabs
  int finestLevel = amrDataPtrArray[0]->FinestLevel();
  Box probDomain(amrDataPtrArray[0]->ProbDomain()[finestLevel]);                                  // Declare a box that's the size of the domain // This is the one we need to check
  int probDomainHeight = probDomain.length(ZDIR);
  int halfProbDomainWidth = (probDomain.length(XDIR))>>1;
  Real probDomainXLength = 1.; // Fix me
  Real probDomainYLength = 1.;
  Real probDomainZLength = 1.;
  Real probDomainVolume  = probDomainXLength * probDomainYLength * probDomainZLength;
    
  // Get some sizes and stuff
  const Real *dx = amrDataPtrArray[0]->DxLevel()[finestLevel].dataPtr();
  VSHOWVAL(verbose, *dx);

  int ix = probDomain.length(XDIR);
  int jx = probDomain.length(YDIR);
  int kx = probDomain.length(ZDIR);

#ifdef DO_SPECTRA
  int wavenumbers = evalWavenumbers(ix,jx,kx); // assume cubic for now
  int hx[DIMS];   hx[0] = ix>>1;   hx[1] = jx>>1;   hx[2] = kx>>1;
#ifdef DO_SPECTRA_PARALLEL
  // Let's try to do a parallel fft (Have to use v2.1.5 of FFTW)

  // Call to mpi plan tells us which slab should be on this processor
  //   Not to rock the boat, need to tell amrlib which box we want to read in here

  // Outline:
  // Call plan
  // Set up which box is on this processor
  // Loop through plots
  // Call FFT
  int local_kx, local_k_start;
  int local_jx_after_transpose;
  int local_j_start_after_transpose;
  int total_local_size;

  fftwnd_mpi_plan plan, inverse_plan;
  fftw_complex *local_xdata, *local_ydata, *local_zdata;

  plan         = fftw3d_mpi_create_plan(MPI_COMM_WORLD,
					kx, jx, ix,
					FFTW_FORWARD, FFTW_ESTIMATE);

  inverse_plan = fftw3d_mpi_create_plan(MPI_COMM_WORLD,
					kx, jx, ix,
					FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwnd_mpi_local_sizes(plan,
			 &local_kx,
			 &local_k_start,
			 &local_jx_after_transpose,
			 &local_j_start_after_transpose,
			 &total_local_size);

  // More fft stuff
  local_xdata = (fftw_complex*) malloc(sizeof(fftw_complex) * total_local_size);
  local_ydata = (fftw_complex*) malloc(sizeof(fftw_complex) * total_local_size);
  local_zdata = (fftw_complex*) malloc(sizeof(fftw_complex) * total_local_size);

  // Now sort boxes as prescribed by FFTW_MPI
  Array<int> tempBoxSmall(nProcs,0);
  Array<int> tempBoxBig(nProcs,0);

  //  std::cout << "MyProc: " << myProc << std::endl;

  tempBoxSmall[myProc] = local_k_start;
  tempBoxBig[myProc]   = tempBoxSmall[myProc] + local_kx - 1;

  for (int iProc=0; iProc<nProcs; iProc++) {
    ParallelDescriptor::ReduceIntSum(tempBoxSmall[iProc]);
    ParallelDescriptor::ReduceIntSum(tempBoxBig[iProc]);
  }
    
  // Divide the domain into slabs as prescribed by FFTW_MPI
  BoxArray domainBoxArray(nProcs);                                                                // Declare an array of boxes to hold our slabs
  Box tempBox(probDomain);                                                                        // Make a box the size of the domain

  for (int iProc=0; iProc<nProcs; iProc++) {
    tempBox.setSmall(ZDIR, tempBoxSmall[iProc]);                                                  // Set new box lo z
    tempBox.setBig(ZDIR, tempBoxBig[iProc]);                                                      // Set new box hi z
    domainBoxArray.set(iProc, tempBox);                                                           // Set size of iProc-th box equal to tempBox
  }
#else
  fftw_plan plan, inverse_plan;
  fftw_complex *local_xdata, *local_ydata, *local_zdata;

  local_xdata = (fftw_complex*) malloc(sizeof(fftw_complex) * ix*jx*kx);
  local_ydata = (fftw_complex*) malloc(sizeof(fftw_complex) * ix*jx*kx);
  local_zdata = (fftw_complex*) malloc(sizeof(fftw_complex) * ix*jx*kx);
  /*
  plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_FORWARD, FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
  plan = fftw_plan_dft_3d(kx,jx,ix,local_xdata,local_xdata,FFTW_BACKWARD,FFTW_MEASURE);   fftw_execute(plan);    fftw_destroy_plan(plan);
  */
  // Make boxes without FFTW
  int iProc=0;
  BoxArray domainBoxArray(nProcs);                                                                // Declare an array of boxes to hold our slabs
  Box tempBox(probDomain);                                                                        // Make a box the size of the domain
  domainBoxArray.set(iProc, tempBox);                                                             // Set size of iProc-th box equal to tempBox
#endif
#endif
  ;                                                                                               // Now let's make sure that this processor gets the right box
  Array<int> pmap(nProcs+1);
  for (int iProc=0; iProc<nProcs; iProc++)
    pmap[iProc] = iProc;
  pmap[nProcs] = myProc;
  DistributionMapping domainDistMap(pmap);                                                        // Declare a distribution mapping

  // Initialise the variables we want to work with, and destination components                    // Would be nice to do this from an inputs file // How do we handle derived variables?
#ifdef DO_FORCING
  if (ParallelDescriptor::IOProcessor())
    std::cout << "DOING FORCING" << std::endl;
  const int nVars(6);
#else
  const int nVars(3);
#endif
  Array<string> whichVar(nVars);
  Array<int>    destFills(nVars);
#ifdef MAESTRO
  whichVar[0] = "x_vel";
  whichVar[1] = "y_vel";
  whichVar[2] = "z_vel";
#else
  whichVar[0] = "x_velocity";
  whichVar[1] = "y_velocity";
  whichVar[2] = "z_velocity";
#endif
#ifdef DO_FORCING
  whichVar[3] = "forcex";
  whichVar[4] = "forcey";
  whichVar[5] = "forcez";
#endif
  for (int c=0; c<nVars; c++ ) destFills[c] = c;                                                  // Destination components required for FillVar()

  // Create a MultiFab with our new boxes, to hold nVars components
  int ngrow(1);
  MultiFab mf;
  mf.define(domainBoxArray, nVars, ngrow, domainDistMap, Fab_allocate);

  // Now let's do all the bits for filling the ghost cells
  int coordSys(0);
  Array<Real> prob_lo = amrDataPtrArray[0]->ProbLo();
  Array<Real> prob_hi = amrDataPtrArray[0]->ProbHi();
  Array<Box> prob_domain = amrDataPtrArray[0]->ProbDomain();
  RealBox rb(prob_lo.dataPtr(),prob_hi.dataPtr());
  Geometry geom;
  geom.define(prob_domain[finestLevel], &rb, coordSys);

  // Initialise the data we want to derive
  Array<Real> AverageData(nPlotFiles*BARVARS, 0.0);                                               // Make an array for each average variable at each time
#ifdef DO_SPECTRA
  Array<Real> SpectraData(nPlotFiles*wavenumbers*DIMS, 0.0);                                      // Make an array for each spectrum at each time
  Real ****CorrltnData = mallocCorrltnData(nPlotFiles,ix>>1,jx>>1,kx>>1);
  Real **IntLen = (Real **) malloc(nPlotFiles*sizeof(Real*));
  Real **Klmgrv = (Real **) malloc(nPlotFiles*sizeof(Real*));
  for (int i=0; i<nPlotFiles; i++) {
      IntLen[i] = (Real*) malloc(DIMS*sizeof(Real));
      Klmgrv[i] = (Real*) malloc(DIMS*sizeof(Real));
  }
#endif

  int **PdfData = (int**) malloc(nPlotFiles*sizeof(int*));
  for (int iPlot = 0; iPlot<nPlotFiles; iPlot++) {
      PdfData[iPlot] = (int*)  malloc(12*PDFBINS*sizeof(int));
      for (int bin=0; bin<12*PDFBINS; bin++)
	  PdfData[iPlot][bin] = 0;
  }
  Real **PdfX = (Real**) malloc(nPlotFiles*sizeof(Real*));
  for (int iPlot = 0; iPlot<nPlotFiles; iPlot++)
      PdfX[iPlot] = (Real*)  malloc(12*PDFBINS*sizeof(Real));
  Real PdfMin[12], PdfMax[12];
  
  // Declare some space for filenames
  char *filename     = (char*) malloc(128*sizeof(char));

  // Let's rock...
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {
    if (verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << std::endl << "--" << std::endl;
      std::cout << "iPlot " << iPlot << "/" << nPlotFiles << std::endl;
    }

    if (verbose && ParallelDescriptor::IOProcessor()) std::cout << " + FillVar()" << std::endl
								<< "   + finestLevel " << finestLevel << std::endl
								<< "   + whichVar " << whichVar[0] << " " << whichVar[1] << " " << whichVar[2] 
#ifdef DO_FORCING
								<< " " << whichVar[3] << " " << whichVar[4] << " " << whichVar[5] 
#endif
								<< " " << std::endl
								<< "   + destFills " <<  destFills[0] << " " << destFills[1] << " " << destFills[2] 
#ifdef DO_FORCING
								<< " " << destFills[3] << " " << destFills[4] << " " << destFills[5] 
#endif
								<< " " << std::endl;
    amrDataPtrArray[iPlot]->FillVar(mf, finestLevel, whichVar, destFills);                        // Populate the MultiFab  

    for (int n=0; n<nVars; n++)
      amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(whichVar[n]));       // Flush grids

    // Get ghost cells if required
    if (ngrow) {
      for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
	const Box& box = mfi.validbox();
	int ncomp = mf.nComp();
	FORT_PUSHVTOGFO(box.loVect(), box.hiVect(),
			mf[mfi].dataPtr(),
			ARLIM(mf[mfi].loVect()),
			ARLIM(mf[mfi].hiVect()),&ncomp);
      }
      mf.FillBoundary(geom.periodicity());
    }

#ifdef MODIFY_SPECTRA
    if (verbose && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calling Spectra... " << std::endl;
    ModifySpectra(*(amrDataPtrArray[iPlot]), mf, plan, inverse_plan,
		  local_xdata, local_ydata, local_zdata);
#endif

    // Get Min/Max for each field
    Real UxMin, UxMax, UyMin, UyMax, UzMin, UzMax, MVMin, MVMax;
#ifdef MAESTRO
    const string xstr = "x_vel";
    const string ystr = "y_vel";
    const string zstr = "z_vel";
    const string vstr = "vort";
#else
    const string xstr = "x_velocity";
    const string ystr = "y_velocity";
    const string zstr = "z_velocity";
    const string vstr = "mag_vort";
#endif
    amrDataPtrArray[iPlot]->MinMax(amrDataPtrArray[iPlot]->ProbDomain()[finestLevel], xstr, 0, UxMin, UxMax);
    amrDataPtrArray[iPlot]->MinMax(amrDataPtrArray[iPlot]->ProbDomain()[finestLevel], ystr, 0, UyMin, UyMax);
    amrDataPtrArray[iPlot]->MinMax(amrDataPtrArray[iPlot]->ProbDomain()[finestLevel], zstr, 0, UzMin, UzMax);
    amrDataPtrArray[iPlot]->MinMax(amrDataPtrArray[iPlot]->ProbDomain()[finestLevel], vstr, 0, MVMin, MVMax);
    PdfMin[0]  = UxMin;                         PdfMax[0]  = UxMax;
    PdfMin[1]  = UyMin;                         PdfMax[1]  = UyMax;
    PdfMin[2]  = UzMin;                         PdfMax[2]  = UzMax;
    MVMax/=2.;
    for (int j=3; j<12; j++) {
	PdfMin[j]  = -MVMax;     PdfMax[j]  = -PdfMin[3];
    }
    for (int j=0; j<12; j++)
	for (int i=0; i<PDFBINS; i++)
	    PdfX[iPlot][j*PDFBINS+i] = PdfMin[j] + (0.5+(Real)i)*(PdfMax[j]-PdfMin[j])/(Real)PDFBINS;
#if 0
    std::cout << "PdfMinMax Ux: " << PdfMin[0] << " " << PdfMax[0] << std::endl;
    std::cout << "PdfMinMax Uy: " << PdfMin[1] << " " << PdfMax[1] << std::endl;
    std::cout << "PdfMinMax Uz: " << PdfMin[2] << " " << PdfMax[2] << std::endl;
    std::cout << "PdfMinMax MV: " << MVMin     << " " << MVMax     << std::endl;
    std::cout << "PdfMinMax DV: " << (UxMax-UxMin)/(2.*dx[0])      << std::endl;
#endif 

    // Call main fortran routine
    if (verbose && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calling Hit... " << std::endl;
    Hit(*(amrDataPtrArray[iPlot]), mf,
	&(AverageData.dataPtr()[iPlot*BARVARS]),
	PdfData[iPlot],PdfMin,PdfMax);                                                            // Process MultiFab data

#ifdef DO_SPECTRA
    if (verbose && ParallelDescriptor::IOProcessor()) 
	std::cout << "Calling Spectra... " << std::endl;
    Spectra(*(amrDataPtrArray[iPlot]), mf, plan, inverse_plan,
	    local_xdata, local_ydata, local_zdata,
	    &(SpectraData.dataPtr()[iPlot*wavenumbers*DIMS]),
	    CorrltnData[iPlot], &(AverageData[iPlot*BARVARS+UX2BAR]));

    if (ParallelDescriptor::IOProcessor()) {

	EvalIntLen(IntLen[iPlot],CorrltnData[iPlot],hx,dx);
	EvalKlmgrv(Klmgrv[iPlot],&(SpectraData.dataPtr()[iPlot*wavenumbers*DIMS]),hx,dx);

	// Write file and check integrals
	Real ux2 = 0.5*AverageData[iPlot*BARVARS+UX2BAR]*probDomainVolume;                          // Get KE from rms
	Real uy2 = 0.5*AverageData[iPlot*BARVARS+UY2BAR]*probDomainVolume;
	Real uz2 = 0.5*AverageData[iPlot*BARVARS+UZ2BAR]*probDomainVolume;

	Real Ex, Ey, Ez;
	Ex = Ey = Ez = 0.;

	sprintf(filename,"%s/SpectralData.dat",plotFileNames[iPlot].c_str());
	std::cout << "Opening " << filename << "..." << std::endl;
	FILE *file = fopen(filename,"w");                                                                 // Output Spectral data
	if (file==NULL)
	  BoxLib::Error("Couldn't open SpectralData.dat");

	for (int wn=0; wn<wavenumbers; wn++) {
	    int i = (iPlot*wavenumbers+wn)*DIMS;
	    fprintf(file,"%e %e %e %e %e %e %e \n",
		    (Real)wn,
		    SpectraData[i],
		    SpectraData[i+1],
		    SpectraData[i+2],
		    pow((Real)wn,(5./3.)) * SpectraData[i],
		    pow((Real)wn,(5./3.)) * SpectraData[i+1],
		    pow((Real)wn,(5./3.)) * SpectraData[i+2]);
	    Ex += SpectraData[i];
	    Ey += SpectraData[i+1];
	    Ez += SpectraData[i+2];
	}
	fclose(file);
	printf("Energy from integration: %e %e %e\n",ux2,uy2,uz2);
	printf("Energy from spectra    : %e %e %e\n",Ex,Ey,Ez);

#ifndef MAESTRO
	sprintf(filename,"%s/QXI.dat",plotFileNames[iPlot].c_str());
	file = fopen(filename,"w");                                                                 // Output x-correlation data
	for (int i=0; i < ix>>1; i++) {
	    fprintf(file,"%e %e %e %e\n", dx[0]*(Real)i, CorrltnData[iPlot][0][0][i], CorrltnData[iPlot][1][0][i], CorrltnData[iPlot][2][0][i]);
	}
	fclose(file);

	sprintf(filename,"%s/QYI.dat",plotFileNames[iPlot].c_str());
	file = fopen(filename,"w");                                                                 // Output y-correlation data
	for (int j=0; j < jx>>1; j++) {
	    fprintf(file,"%e %e %e %e\n", dx[1]*(Real)j, CorrltnData[iPlot][0][1][j], CorrltnData[iPlot][1][1][j], CorrltnData[iPlot][2][1][j]);
	}
	fclose(file);

	sprintf(filename,"%s/QZI.dat",plotFileNames[iPlot].c_str());
	file = fopen(filename,"w");                                                                 // Output z-correlation data
	for (int k=0; k < kx>>1; k++) {
	    fprintf(file,"%e %e %e %e\n", dx[2]*(Real)k, CorrltnData[iPlot][0][2][k], CorrltnData[iPlot][1][2][k], CorrltnData[iPlot][2][2][k]);
	}
	fclose(file);
#endif
    }
#endif

    // Write running files
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Writing data to file..." << std::endl;

	FILE *file;
	sprintf(filename,"%s/AverageData.dat",plotFileNames[iPlot].c_str());
	file = fopen(filename,"w");
	write_average_data(file,&(AverageData[iPlot*BARVARS]),IntLen[iPlot],Klmgrv[iPlot],time[iPlot],dx[0]);
	fclose(file);

#ifndef MAESTRO
	sprintf(filename,"%s/PdfData.dat",plotFileNames[iPlot].c_str());
	file = fopen(filename,"w");                                                                 // Output x-correlation data
	for (int i=0; i < PDFBINS; i++) {
	    for (int v=0; v<12; v++) {
		fprintf(file,"%e %i ", PdfX[iPlot][v*PDFBINS+i], PdfData[iPlot][v*PDFBINS+i]);
	    }
	    fprintf(file,"\n");
	}
	fclose(file);
#endif
    }    
  } // iPlot

#ifdef DO_SPECTRA_PARALLEL
  fftwnd_mpi_destroy_plan(plan);
  fftwnd_mpi_destroy_plan(inverse_plan);
#endif

  if (ParallelDescriptor::IOProcessor()) {                                                        // Finish processing on IOProcessor
    std::cout << "Writing data to file..." << std::endl;
    FILE *file;                                                                                   // Output Average data
    file = fopen("AverageData.dat","w");
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
	write_average_data(file,&(AverageData[iPlot*BARVARS]),IntLen[iPlot],Klmgrv[iPlot],time[iPlot],dx[0]);
    fclose(file);

    file = fopen("AverageData.key","w");
    int i = 1;
    fprintf(file,"%02d time\n",i++);
    fprintf(file,"%02d ux\n",i++);
    fprintf(file,"%02d uy\n",i++);
    fprintf(file,"%02d uz\n",i++);
    fprintf(file,"%02d ux2\n",i++);
    fprintf(file,"%02d uy2\n",i++);
    fprintf(file,"%02d uz2\n",i++);
    fprintf(file,"%02d ux2-ux*ux\n",i++);
    fprintf(file,"%02d uy2-uy*uy\n",i++);
    fprintf(file,"%02d uz2-uz*uz\n",i++);
    fprintf(file,"%02d energy\n",i++);
    fprintf(file,"%02d dissipation (strain rate)\n",i++);
    fprintf(file,"%02d dissipation (Laplacian)\n",i++);
    fprintf(file,"%02d dE/dt\n",i++);
#ifdef DO_FORCING
    fprintf(file,"%02d forcing\n",i++);
    fprintf(file,"%02d beta\n",i++);
#else
    fprintf(file,"%02d forcing (zero)\n",i++);
    fprintf(file,"%02d beta    (zero)\n",i++);
#endif
#ifdef DO_SPECTRA
    fprintf(file,"%02d lx\n",i++);
    fprintf(file,"%02d ly\n",i++);
    fprintf(file,"%02d lz\n",i++);
    fprintf(file,"%02d etax\n",i++);
    fprintf(file,"%02d etay\n",i++);
    fprintf(file,"%02d etaz\n",i++);
#endif
    fclose(file);
    
    std::cout << "   ...done." << std::endl;
  } // IOProcessor for writing AverageData.dat

#ifdef DO_SPECTRA
  freeCorrltnData(CorrltnData,nPlotFiles);
#endif
  BoxLib::Finalize();
  return 0;
    
}



void write_average_data(FILE *file, Real *AverageData, Real *IntLen, Real *Klmgrv, Real time, const Real dx)
{
    Real ux, uy, uz, ux2, uy2, uz2, enrg, diss, dis2, forc, beta;

    ux  = AverageData[UXBAR];    uy  = AverageData[UYBAR];    uz  = AverageData[UZBAR];
    ux2 = AverageData[UX2BAR];   uy2 = AverageData[UY2BAR];   uz2 = AverageData[UZ2BAR];
    enrg= 0.5*(ux2+uy2+uz2);     diss= AverageData[DBAR];     dis2= AverageData[D2BAR];
        
#ifdef DO_FORCING
    forc= AverageData[FBAR];     beta= pow(forc,(2./3.)) / ( diss * pow(dx,(4./3.)) );
#else
    forc= 0.;    beta= 0.;       // Let's leave it in so that there's an entry in the data file
#endif
    //             t ux uy uz x2 y2 z2
    fprintf(file,"%e %e %e %e %e %e %e ",time,ux,uy,uz,ux2,uy2,uz2);
    fprintf(file,"%e %e %e %e %e %e %e %e ", ux2-ux*ux, uy2-uy*uy, uz2-uz*uz, enrg, diss, dis2, forc, beta);
#ifdef DO_SPECTRA
    fprintf(file,"%e %e %e ", IntLen[0], IntLen[1], IntLen[2]);
    fprintf(file,"%e %e %e ", Klmgrv[0], Klmgrv[1], Klmgrv[2]);
#endif
    fprintf(file,"\n");

    return;
}




static
void
PrintUsage (char* progName)
{
  std::cout << "\nUsage:\n"
	    << progName
	    << "\n\tinfile = inputFileName"
	    << "\n\t[-help]"
	    << "\n\n";
  exit(1);
}






void qs_swap(int *a, int *b) { int t=*a; *a=*b; *b=t; }

void qs_sort(int arr[], int beg, int end)
{
  if (end > beg + 1) {
    int piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else
	qs_swap(&arr[l], &arr[--r]);
    }
    qs_swap(&arr[--l], &arr[beg]);
    qs_sort(arr, beg, l);
    qs_sort(arr, r, end);
  }
}




Real ****mallocCorrltnData(int nPlotFiles, int ix, int jx, int kx)
{
  Real ****CorrltnData;

  CorrltnData = (Real****)malloc(nPlotFiles*sizeof(Real***));
  if (CorrltnData==NULL) std::cout << "Argh..." << std::endl;

  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {
    CorrltnData[iPlot] = (Real***)malloc(3*sizeof(Real**));     
    if (CorrltnData[iPlot]==NULL) std::cout << "Argh..." << std::endl;
    
    for (int c=0; c<3; c++) {
      CorrltnData[iPlot][c] = (Real**)malloc(3*sizeof(Real*));  
      if (CorrltnData[iPlot][c]==NULL) std::cout << "Argh..." << std::endl;

      CorrltnData[iPlot][c][0] = (Real*)malloc(ix*sizeof(Real)); for (int i=0; i<ix; i++) CorrltnData[iPlot][c][0][i] = 0.;
      if (CorrltnData[iPlot][c][0]==NULL) std::cout << "Argh..." << std::endl;

      CorrltnData[iPlot][c][1] = (Real*)malloc(jx*sizeof(Real)); for (int j=0; j<jx; j++) CorrltnData[iPlot][c][1][j] = 0.;
      if (CorrltnData[iPlot][c][1]==NULL) std::cout << "Argh..." << std::endl;

      CorrltnData[iPlot][c][2] = (Real*)malloc(kx*sizeof(Real)); for (int k=0; k<kx; k++) CorrltnData[iPlot][c][2][k] = 0.;
      if (CorrltnData[iPlot][c][2]==NULL) std::cout << "Argh..." << std::endl;
    }
  }

  return(CorrltnData);
}





void freeCorrltnData(Real ****CorrltnData, int nPlotFiles)
{
  for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	free(CorrltnData[iPlot][i][j]);
      }
      free(CorrltnData[iPlot][i]);
    }
    free(CorrltnData[iPlot]);
  }
  free(CorrltnData);
}


int evalWavenumbers(int ix, int jx, int kx)
{
    int wn;

    if (ix>jx)
	wn = ix;
    else
	wn = jx;
    if (kx>wn)
	wn = kx;

    return(wn>>1);
}







#if 0
    int inverse_local_kx, inverse_local_k_start;
    int inverse_local_jx_after_transpose;
    int inverse_local_j_start_after_transpose;
    int inverse_total_local_size;

    fftwnd_mpi_local_sizes(inverse_plan,
			   &inverse_local_kx,
			   &inverse_local_k_start,
			   &inverse_local_jx_after_transpose,
			   &inverse_local_j_start_after_transpose,
			   &inverse_total_local_size);

    int fore, back;
    fore = local_kx;  back = inverse_local_kx;
    if (fore != back) std::cout << ParallelDescriptor::MyProc() << ": local_kx (" << fore << ", " << back << ")." << std::endl;
    fore = local_k_start;  back = inverse_local_k_start;
    if (fore != back) std::cout << ParallelDescriptor::MyProc() << ": local_k_start (" << fore << ", " << back << ")." << std::endl;
    fore = local_jx_after_transpose;  back = inverse_local_jx_after_transpose;
    if (fore != back) std::cout << ParallelDescriptor::MyProc() << ": local_jx_after_transpose (" << fore << ", " << back << ")." << std::endl;
    fore = local_j_start_after_transpose;  back = inverse_local_j_start_after_transpose;
    if (fore != back) std::cout << ParallelDescriptor::MyProc() << ": local_j_start_after_transpose (" << fore << ", " << back << ")." << std::endl;
    fore = total_local_size;  back = inverse_total_local_size;
    if (fore != back) std::cout << ParallelDescriptor::MyProc() << ": total_local_size (" << fore << ", " << back << ")." << std::endl;
#endif







#if 0

    for (int k = 0; k < kx; k++) {
      for (int j = 0; j < jx; j++) {
	for (int i = 0; i < ix; i++) {
	  int cell = (k*jx + j)*ix + i;
	  local_xdata[cell].re = ux[cell];   local_xdata[cell].im = 0.;
	  local_ydata[cell].re = uy[cell];   local_ydata[cell].im = 0.;
	  local_zdata[cell].re = uz[cell];   local_zdata[cell].im = 0.;
	}
      }
    }
    fftwnd_mpi(plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);

    fftwnd_mpi(inverse_plan, 1, local_xdata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_ydata, NULL, FFTW_NORMAL_ORDER);
    fftwnd_mpi(inverse_plan, 1, local_zdata, NULL, FFTW_NORMAL_ORDER);

    int cell;
    Real ur, ui, vr, vi, wr, wi;
    Real xsum = 0.;
    Real ysum = 0.;
    Real zsum = 0.;
    for (int k = 0; k < kx; k++) {
      for (int j = 0; j < jx; j++) {
	for (int i = 0; i < ix; i++) {
	  cell = (k*jx + j)*ix + i;
	  ur = local_xdata[cell].re / (Real)(dix*djx*dkx) - ux[cell];
	  ui = local_xdata[cell].im / (Real)(dix*djx*dkx);
	  vr = local_ydata[cell].re / (Real)(dix*djx*dkx) - uy[cell];
	  vi = local_ydata[cell].im / (Real)(dix*djx*dkx);
	  wr = local_zdata[cell].re / (Real)(dix*djx*dkx) - uz[cell];
	  wi = local_zdata[cell].im / (Real)(dix*djx*dkx);
	  xsum += ur*ur + ui*ui;
	  ysum += vr*vr + vi*vi;
	  zsum += wr*wr + wi*wi;
	}
      }
      std::cout << "ur " << ur << std::endl;
      std::cout << "ui " << ui << std::endl;
      std::cout << "vr " << vr << std::endl;
      std::cout << "vi " << vi << std::endl;
      std::cout << "wr " << wr << std::endl;
      std::cout << "wi " << wi << std::endl;
    }
    ParallelDescriptor::ReduceRealSum(xsum);
    ParallelDescriptor::ReduceRealSum(ysum);
    ParallelDescriptor::ReduceRealSum(zsum);
    VSHOWVAL(1,xsum);
    VSHOWVAL(1,ysum);
    VSHOWVAL(1,zsum);

#endif
