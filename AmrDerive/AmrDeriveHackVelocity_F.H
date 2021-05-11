#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_HACKVELOCITY          hackvelocity
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_HACKVELOCITY          HACKVELOCITY
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_HACKVELOCITY          hackvelocity
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_HACKVELOCITY          hackvelocity_
#  endif

//
// Use ARLIM in the following.
//

extern "C"
{

  void FORT_HACKVELOCITY (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		    Real *outData, ARLIM_P(outDat_lo), ARLIM_P(outDat_hi),
		    const Real* probLo, const Real* probHi, const Real* delta, 
		    const int *nFiltVars, const int *fSize);

} 

#endif
