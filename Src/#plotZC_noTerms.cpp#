#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

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
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
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
      AmrData::SetVerbose(false);

    std::string plotFileName; pp.get("infile",plotFileName);
    std::string fuelName = "H2"; pp.query("fuelName",fuelName);
    Real s = 8.0; pp.query("stoichRatio",s);
    Real Y_O_air = 0.233; pp.query("YO2Air",Y_O_air);
    Real S = s/Y_O_air;
    const Real Z_st = Y_O_air/(s+Y_O_air);
    Real Tmin = 300; pp.query("Tmin",Tmin);
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    //init_mech();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    
    int idFCRin = -1;
    int idFin = -1;
    int idOin = -1;
    int idRhoin = -1;
    int idTin = -1;
    int idVortin = -1;
    //Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string FCRName = "I_R("+fuelName+")";//"_ConsumptionRate";
    const std::string spName= "Y("+fuelName+")";
    const std::string oxName= "Y(O2)";
    const std::string rhoName = "density";
    const std::string tName = "temp";
    const std::string vortName = "mag_vort";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == FCRName) idFCRin = i;
      if (plotVarNames[i] == spName) idFin = i;
      if (plotVarNames[i] == oxName) idOin = i;
      if (plotVarNames[i] == rhoName) idRhoin = i;
      if (plotVarNames[i] == tName) idTin = i;
      if (plotVarNames[i] == vortName) idVortin = i;
    }
    if (idFin<0 || idOin<0 ||  idRhoin < 0 || idTin < 0 || idFCRin < 0 || idVortin < 0)
      Abort("Cannot find required data in pltfile");
    const int nCompIn  = 6;
    int nCompOut = 4;
    
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    //in
    const int idFlocal = 0; // fuel in here
    const int idOlocal = 1; // O2 in here
    const int idRholocal = 2; // density in here
    const int idFCRlocal = 3; //FCR in here
    const int idTlocal = 4; //Temp in here
    const int idVortlocal = 5;
    //out
    const int idZlocal = 0; // Z out here
    const int idClocal = 1; // C out here
    const int idCsourcelocal = 2; //C source here
    const int idDalocal = 3;
    //const int idPhilocal = 4;
    destFillComps[idFlocal] = idFlocal;
    destFillComps[idOlocal] = idOlocal;
    destFillComps[idRholocal] = idRholocal;
    destFillComps[idFCRlocal] = idFCRlocal;
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idVortlocal] = idVortlocal;
    inNames[idFlocal] = spName;
    inNames[idOlocal] =  oxName;
    inNames[idRholocal] = rhoName;
    inNames[idFCRlocal] = FCRName;
    inNames[idTlocal] = tName;
    inNames[idVortlocal] = vortName;
    outNames[idZlocal] = "mixture_fraction";
    outNames[idClocal] = "progress_variable";
    outNames[idCsourcelocal] = "progress_variable_source";
    outNames[idDalocal] = "log(Da)";
    //outNames[idPhilocal] = "phi";
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<MultiFab*> baseComps(Nlev);
    Vector<MultiFab*> indata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    const int nGrow = 1;
    int b[3] = {1, 1, 1};
    Real eps = 1e-8;
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const Vector<Real>& delta = amrData.DxLevel()[lev];
      const DistributionMapping dm(ba);
      indata[lev] = new MultiFab(ba,dm,nCompIn,nGrow);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(*indata[lev],lev,inNames,destFillComps); //Problem
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,0,&(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
        Array4<Real> const& ZC = (*outdata[lev]).array(mfi);
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  
	  ZC(i,j,k,idZlocal) = (s*inbox(i,j,k,idFlocal) - inbox(i,j,k,idOlocal) + Y_O_air)/(s+Y_O_air);
	  if (inbox(i,j,k,idTlocal) < Tmin || ZC(i,j,k,idZlocal) < eps) {
	    ZC(i,j,k,idClocal) = 0;
	  } else if (ZC(i,j,k,idZlocal) <= Z_st) {
	    ZC(i,j,k,idClocal) = 1-inbox(i,j,k,idFlocal)/ZC(i,j,k,idZlocal);
	    ZC(i,j,k,idCsourcelocal) = -inbox(i,j,k,idFCRlocal)/ZC(i,j,k,idZlocal);
	  } else {
	    ZC(i,j,k,idClocal) = ((1-Z_st)/Z_st)*(ZC(i,j,k,idZlocal)-inbox(i,j,k,idFlocal))/(1-ZC(i,j,k,idZlocal));
	    ZC(i,j,k,idCsourcelocal) = -((1-Z_st)/Z_st)*inbox(i,j,k,idFCRlocal)/(1-ZC(i,j,k,idZlocal));
	  }
	  ZC(i,j,k,idDalocal) = std::log10(std::max(ZC(i,j,k,idCsourcelocal)/(inbox(i,j,k,idRholocal)*inbox(i,j,k,idVortlocal)),1e-10));
	  //ZC(i,j,k,idPhilocal) = S*ZC(i,j,k,idZlocal)/(1-ZC(i,j,k,idZlocal));
        });
	
      }
      
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_ZC");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
