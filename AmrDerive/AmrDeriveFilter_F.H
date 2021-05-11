#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_FILTER          filter
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_FILTER          FILTER
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_FILTER          filter
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_FILTER          filter_
#  endif

//
// Use ARLIM in the following.
//

extern "C"
{

  void FORT_FILTER (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		    Real *outData, ARLIM_P(outDat_lo), ARLIM_P(outDat_hi),
		    const Real* probLo, const Real* probHi, const Real* delta, 
		    const int *nFiltVars, const int *fSize, 
		    const int *favre, const int *denVar);

} 

#endif
