#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>
#include <mechanism.h>
#include <chemistry_file.H>
#include <util.H>
#include <util_F.H>
#include <Transport_F.H>
#include <Fuego_EOS.H>

using namespace amrex;
using namespace analysis_util;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    init_mech();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    int idRin = -1;
    int idTin = -1;
    int idQin = -1; 
    Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string RName = "density";
    const std::string TName = "temp";
    const std::string QName = "HeatRelease";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == RName)  idRin = i;
      if (plotVarNames[i] == TName)  idTin = i;
      if (plotVarNames[i] == QName)  idQin = i;
    }
    if (idRin<0 || idYin<0 || idTin<0 || idQin < 0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int nCompIn  = NUM_SPECIES + 3;
    const int idLeout   = 0;
    const int nCompOut = 1;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idRlocal = NUM_SPECIES;   // Rho start here
    const int idTlocal = NUM_SPECIES+1; // T starts here
    const int idQlocal = NUM_SPECIES+2; // Q starts here
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    destFillComps[idRlocal] = idRlocal;
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idQlocal] = idQlocal;
    inNames[idRlocal] = RName;
    inNames[idTlocal] = TName;
    inNames[idQlocal] = QName;
    outNames[0] = "QoverRhoCP";

    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;
    int b[3] = {1, 1, 1};

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& Y  = indata.array(mfi);
        Array4<Real> const& R  = indata.array(mfi);
        Array4<Real> const& T  = indata.array(mfi);
        Array4<Real> const& Q  = indata.array(mfi);
        Array4<Real> const& Qcp = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
          Real Yl[NUM_SPECIES];
          for (int n=0; n<NUM_SPECIES; ++n) {
            Yl[n] = Y(i,j,k,idYlocal+n);
          }
          Real Cpmix;
          CKCPBS(&T(i,j,k,idTlocal),Yl,&Cpmix);
	  Cpmix*=1.e-4;
	  Qcp(i,j,k,0) = Q(i,j,k,idQlocal)/R(i,j,k,idRlocal)/Cpmix;
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_QoverRhoCP");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
