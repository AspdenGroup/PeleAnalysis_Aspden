// --------------------------------------------------------------------
// CoarsenPLT.cpp
// --------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <DataServices.H>
#include <Utility.H>
#include <FArrayBox.H>
#include <Utility.H>
#include <AmrDeriveCoarsen_F.H>
#include <WritePlotFile.H>


// --------------------------------------------------------------------
static void print_usage(int, char *argv[]) {
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<filename> [numNewCoarseLevels=<n> [DEF=1] ";
  std::cerr << "finestLevel=<n> [DEF=infile.finestLev]]\n";
    exit(1);
}


// --------------------------------------------------------------------
void avgDown_doit(const FArrayBox &fine_fab,
             FArrayBox &crse_fab, 
             const Box &ovlp,
             int        scomp,
             int        dcomp,
             int        ncomp,
             int        ratio)
{
    const int  *ovlo   = ovlp.loVect();
    const int  *ovhi   = ovlp.hiVect();
    const int  *flo    = fine_fab.loVect();
    const int  *fhi    = fine_fab.hiVect();
    const Real *f_dat  = fine_fab.dataPtr(scomp);
    const int  *clo    = crse_fab.loVect();
    const int  *chi    = crse_fab.hiVect();
    Real       *c_dat  = crse_fab.dataPtr(dcomp);

    FORT_AVGDN(ovlo,ovhi,
               c_dat,ARLIM(clo),ARLIM(chi),
               f_dat,ARLIM(flo),ARLIM(fhi),
               &ratio,&ncomp);
}


// --------------------------------------------------------------------
void avgDown(MultiFab  &S_crse,
        MultiFab       &S_fine,
        int             scomp,
        int             dcomp,
        int             ncomp,
        int             ratio)
{
    const BoxArray &fgrids = S_fine.boxArray();
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA(fgrids.size());

    for (int i = 0; i < fgrids.size(); ++i)
    {
        crse_S_fine_BA.set(i,BoxLib::coarsen(fgrids[i],ratio));
    }

    MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);

    for(MFIter mfi(S_fine); mfi.isValid(); ++mfi) {
        const int i(mfi.index());
        avgDown_doit(S_fine[i],crse_S_fine[i],
                     crse_S_fine_BA[i],scomp,0,ncomp,ratio);
    }
    S_crse.copy(crse_S_fine,0,dcomp,ncomp);
}


// --------------------------------------------------------------------
int main(int argc, char *argv[]) {
    BoxLib::Initialize(argc,argv);

    if(argc < 2) {
      print_usage(argc,argv);
    }

    ParmParse pp;

    if(pp.contains("help")) {
      print_usage(argc,argv);
    }

    int verbose = 0; pp.query("verbose",verbose);
    if(verbose > 1) {
      AmrData::SetVerbose(true);
    }

    std::string infile; pp.get("infile",infile);
    std::string outfile(infile + "_crse"); pp.query("outfile",outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData &amrData = dataServices.AmrDataRef();

    int finestLevel(amrData.FinestLevel());
    pp.query("finestLevel",finestLevel);

    int numNewCoarseLevels = 1; pp.query("numNewCoarseLevels",numNewCoarseLevels);

    Array<int> comps;
    if(int nc = pp.countval("comps")) {
      comps.resize(nc);
      pp.getarr("comps",comps,0,nc);
    } else {
      int sComp = 0; pp.query("sComp",sComp);
      int nComp = amrData.PlotVarNames().size(); pp.query("nComp",nComp);
      BL_ASSERT(sComp+nComp <= amrData.NComp());
      comps.resize(nComp);
      for(int i(0); i < nComp; ++i) {
        comps[i] = sComp + i;
      }
    }

    int Nlev(finestLevel+1+numNewCoarseLevels);
    PArray<MultiFab> state(Nlev);
    const Array<string>& plotVarNames = amrData.PlotVarNames();
    
    int numState(comps.size());
    Array<string> inVarNames(numState);
    Array<int> destFillComps(numState);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "Filling the following:";
    }
    for(int i(0); i < numState; ++i) {
      destFillComps[i] = comps[i];
      inVarNames[i] = plotVarNames[comps[i]];
      if(ParallelDescriptor::IOProcessor()) {
        cout << " " << inVarNames[i];
      }
    }
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl;
    }
    int nGrow(0);

    Array<int> refRatio(Nlev-1,2);
    Array<Box> probDomain(Nlev);
    Array< Array<Real> > dxLevel(Nlev);
    for(int lev(0); lev <= finestLevel; ++lev) {
      int shiftedLevel = lev + numNewCoarseLevels;
      if(lev < finestLevel) {
        refRatio[shiftedLevel] = amrData.RefRatio()[lev];
      }
      probDomain[shiftedLevel] = amrData.ProbDomain()[lev];
      dxLevel[shiftedLevel].resize(BL_SPACEDIM);
      for(int i(0); i < BL_SPACEDIM; ++i) {
        dxLevel[shiftedLevel][i] = amrData.DxLevel()[lev][i];
      }

      const BoxArray ba = amrData.boxArray(lev);
      state.set(shiftedLevel,new MultiFab(ba,numState,nGrow));

      amrData.FillVar(state[shiftedLevel],lev,inVarNames,destFillComps);
      for(int i(0); i < inVarNames.size(); ++i) {
        amrData.FlushGrids(amrData.StateNumber(inVarNames[i]));
      }

      if(ParallelDescriptor::IOProcessor()) {
       cerr << "Data has been read for level " << lev << endl;
      }
    }
    for(int lev(-1); lev >= -numNewCoarseLevels; --lev) {
      int myLev = numNewCoarseLevels+lev;
      const BoxArray ba = BoxArray(state[myLev+1].boxArray()).coarsen(refRatio[myLev]);
      probDomain[myLev] = Box(probDomain[myLev+1]).coarsen(refRatio[myLev]);
      dxLevel[myLev].resize(BL_SPACEDIM);
      for(int i(0); i < BL_SPACEDIM; ++i) {
        dxLevel[myLev][i] = refRatio[myLev] * dxLevel[myLev+1][i];
      }
      state.set(myLev,new MultiFab(ba,numState,nGrow));
      avgDown(state[myLev],state[myLev+1],0,0,numState,refRatio[myLev]);
    }

    if(ParallelDescriptor::IOProcessor()) {
      cout << "Writing new data to " << outfile << endl;
    }
    
    bool verb(false);
    const AmrData &adat = amrData;
    std::string plotFileVersion = "NavierStokes-V1.1";
    WritePlotfile(plotFileVersion, state, adat.Time(),
                  adat.ProbLo(), adat.ProbHi(), refRatio,
                  probDomain, dxLevel, adat.CoordSys(),
		  outfile, inVarNames, verb);

    BoxLib::Finalize();

    return 0;
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------
