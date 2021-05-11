#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_INTEGRATE          integrate
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_INTEGRATE          INTEGRATE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_INTEGRATE          integrate
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_INTEGRATE          integrate_
#  endif

//
// Use ARLIM in the following.
//

#define DIMS 3

extern "C"
{
  void FORT_INTEGRATE (const Real* dat, int *nVars,
			ARLIM_P(dat_lo), ARLIM_P(dat_hi),
			ARLIM_P(grid_lo), ARLIM_P(grid_hi),
			const Real* delta, Real* turb,
			int *isize, int *jsize, int *ksize, int *turbVars);
} 

#define VSHOWVAL(verbose, val) { if (verbose && ParallelDescriptor::IOProcessor()) { \
                 std::cout << #val << " = " << val << '\n'; } }

static void PrintUsage (char* progName);

#endif
