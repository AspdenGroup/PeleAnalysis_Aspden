#define XDIR 0
#define YDIR 1
#define ZDIR 2
#define DIM  3

int  nProcs;
int  myProc;
int  IOProc;

int  verbose;

int  nVars;

int  div_free;
int  transpose_dp;
int  density_weighting;

int  BLix, BLjx, BLkx;
int  FTix, FTjx, FTkx;
int  FTis, FTjs, FTks;
int  FTmx, FThkxpo;
int  wavenumbers;

int  Qix, Qjx, Qkx;
Array<Real*> spectrum;
Array<Real*> spectrumS;
Array<Real*> spectrumC;
Array<Real*> Qx, Qy, Qz;
Real *sum, *sum2;

Box probDomain;
Array<Real> probLo;
Array<Real> probHi;
Real Lx, Ly, Lz;
Real dx, dy, dz;

std::string infile;
Array<std::string> whichVar;
Real Time;
int  timeSteps;

int  local_ix, local_i_start;
int  local_jx_after_transpose, local_j_start_after_transpose;
int  total_local_size;
Array<Real> local_xlo;
Array<Real> local_xhi;

Array<fftw_real*>    local_data;
Array<fftw_complex*> local_data_c;

rfftwnd_mpi_plan plan_real2cplx, plan_cplx2real; 

void plan_ffts(Box &probDomain, BoxArray &domainBoxArray, Array<int> &pmap);

void Spectra(MultiFab &mf);

int        do_filter;
int        nFilters;
Array<int> filterWN;
void       filter(MultiFab &mf);
