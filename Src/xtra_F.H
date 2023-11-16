
#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_PUSHVTOGHO      pushvtogho
#  define   FORT_PUSHVTOGFO      pushvtogfo
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_PUSHVTOGHO      PUSHVTOGHO
#    define   FORT_PUSHVTOGFO      PUSHVTOGFO
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_PUSHVTOGHO      pushvtogho
#    define   FORT_PUSHVTOGFO      pushvtogfo
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_PUSHVTOGHO      pushvtogho_
#    define   FORT_PUSHVTOGFO      pushvtogfo_
#  endif

//
// Use ARLIM in the following.
//

extern "C"
{
    void FORT_PUSHVTOGHO(const int* dlo, const int* dhi,
                         Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                         const int* nc);

    void FORT_PUSHVTOGFO(const int* dlo, const int* dhi,
                         Real* U, ARLIM_P(Ulo), ARLIM_P(Uho),
                         const int* nc);
}

#endif
