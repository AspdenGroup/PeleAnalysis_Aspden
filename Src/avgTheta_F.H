#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_INTEGRATE          integrate
#  define   FORT_PUSHVTOGHO      pushvtogho
#  define   FORT_PUSHVTOGFO      pushvtogfo
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_INTEGRATE          INTEGRATE
#    define   FORT_PUSHVTOGHO      PUSHVTOGHO
#    define   FORT_PUSHVTOGFO      PUSHVTOGFO
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_INTEGRATE          integrate
#    define   FORT_PUSHVTOGHO      pushvtogho
#    define   FORT_PUSHVTOGFO      pushvtogfo
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_INTEGRATE          integrate_
#    define   FORT_PUSHVTOGHO      pushvtogho_
#    define   FORT_PUSHVTOGFO      pushvtogfo_
#  endif

//
// Use ARLIM in the following.
//

#define DIMS 3

extern "C"
{
  void FORT_INTEGRATE (const Real* dat,
			ARLIM_P(dat_lo), ARLIM_P(dat_hi),
			ARLIM_P(grid_lo), ARLIM_P(grid_hi),
			const Real* xlo, const Real* delta, Real* turb,
			int *rsize, int *levRsize,int *ksize, int *levKsize,
		        int *nVars, int *turbVars);
  void FORT_PUSHVTOGHO(const int* dlo, const int* dhi,
                         Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                         const int* nc);
  void FORT_PUSHVTOGFO(const int* dlo, const int* dhi,
                         Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                         const int* nc);
} 

#define VSHOWVAL(verbose, val) { if (verbose && ParallelDescriptor::IOProcessor()) { \
                 std::cout << #val << " = " << val << '\n'; } }

static void PrintUsage (char* progName);

#endif
