/* -------------------------------------------------------------------- */
/*  CoarsenPLT_F.H */
/* -------------------------------------------------------------------- */
#ifndef _AMRDERIVECOARSEN_F_H_
#define _AMRDERIVECOARSEN_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_AVGDN     avgdn
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_AVGDN     AVGDN
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_AVGDN     avgdn
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_AVGDN     avgdn_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_AVGDN(const int* lo, const int* hi,
                    Real*       c, ARLIM_P(clo), ARLIM_P(chi),
                    const Real* f, ARLIM_P(flo), ARLIM_P(fhi),
                    const int* ratio, const int* nComp);
}
#endif

#endif /*_AMRDERIVECOARSEN */
