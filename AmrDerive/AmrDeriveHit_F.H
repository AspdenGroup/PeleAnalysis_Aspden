#include <REAL.H>
#include <ArrayLim.H>

#define DO_FORCINGundef
#define MAESTRO

#define UXBAR   0
#define UYBAR   1
#define UZBAR   2
#define UX2BAR  3
#define UY2BAR  4
#define UZ2BAR  5
#define DBAR    6
#define D2BAR   7

#ifdef DO_FORCING
#define FBAR    8
#define BARVARS 9
#else
#define BARVARS 8
#endif

#define PDFBINS 1024

#ifdef BL_LANG_FORT
#  define   FORT_HIT          fdhit
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_HIT          FDHIT
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_HIT          fdhit
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_HIT          fdhit_
#  endif

extern "C"
{
  void FORT_HIT (const Real* dat,
		 ARLIM_P(dat_lo),  ARLIM_P(dat_hi),
		 ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		 const Real* delta, Real* avgs,
		 int *PdfData, Real *PdfMin, Real *PdfMax);
} 

#define VSHOWVAL(verbose, val) { if (verbose && ParallelDescriptor::IOProcessor()) { \
                 std::cout << #val << " = " << val << '\n'; } }

#define DIMS 3
#define ITERS 1e4
#define TOL 1.e-8

void write_average_data(FILE *file, Real *AverageData, Real *IntLen, Real *Klmgrv, Real time, const Real dx);

int evalWavenumbers(int ix, int jx, int kx);

static void PrintUsage (char* progName);

static void Hit(AmrData &amrData, const MultiFab &mf, Real *avgs, Real *PdfData, Real *PdfMin, Real *PdfMax);

void qs_swap(int *a, int *b);

void qs_sort(int arr[], int beg, int end);

#endif
