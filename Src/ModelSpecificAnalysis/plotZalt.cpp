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
    Real T0; pp.get("T0",T0);
    Real T1; pp.get("T1",T1);
    Real Y0; pp.get("Y0",Y0);
    Real Y1; pp.get("Y1",Y1);
    int idYin = -1;
    int idTin = -1;
    int idQin = -1;
    int idFCRin = -1;
    Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string TName = "temp";
    const std::string QName = "HeatRelease";
    const std::string FCRName = "I_R(H2)";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName)  idTin = i;
      if (plotVarNames[i] == QName)  idQin = i;
      if (plotVarNames[i] == FCRName) idFCRin = i;
    }
    if (idYin<0 || idTin<0 || idQin < 0 || idFCRin < 0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int nCompIn  = NUM_SPECIES + 3;
    const int nCompOut = 1;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idTlocal = NUM_SPECIES;   // T starts here
    const int idQlocal = NUM_SPECIES+1; // Q starts here
    const int idFCRlocal = NUM_SPECIES+2;
    const int idH2local = -1;
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idQlocal] = idQlocal;
    destFillComps[idFCRlocal] = idFCRlocal;
    inNames[idTlocal] = TName;
    inNames[idQlocal] = QName;
    inNames[idFCRlocal] = FCRName;
    outNames[0] = "Zalt";

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
        Array4<Real> const& T  = indata.array(mfi);
        Array4<Real> const& Q  = indata.array(mfi);
	Array4<Real> const& FCR = indata.array(mfi);
        Array4<Real> const& Z = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
          Real Yl[NUM_SPECIES];
          for (int n=0; n<NUM_SPECIES; ++n) {
            Yl[n] = Y(i,j,k,idYlocal+n);
          }
          Real Cpmix;
	  Real Qscale = -Q(i,j,k,idQlocal)/FCR(i,j,k,idFCRlocal);
          CKCPBS(&T(i,j,k,idTlocal),Yl,&Cpmix);
	  Z(i,j,k,0) = (Cpmix*(T(i,j,k,idTlocal)-T0) + Qscale*(Y(i,j,k,idYlocal+H2_ID)-Y0))/(Cpmix*(T1-T0)+Qscale*(Y1-Y0));
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_Le");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
