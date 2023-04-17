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
    std::string productName = "H2O"; pp.query("productName",productName);
    int clipProgress = 0; pp.query("clipProgress",clipProgress);
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

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    
    int idFPRin = -1;
    int idFin = -1;
    int idOin = -1;
    int idPin = -1;
    int idPPRin = -1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string FPRName = "I_R("+fuelName+")";
    const std::string spName= "Y("+fuelName+")";
    const std::string oxName= "Y(O2)";
    const std::string prodName= "Y("+productName+")";
    const std::string PPRName = "I_R("+productName+")";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == FPRName) idFPRin = i;
      if (plotVarNames[i] == spName) idFin = i;
      if (plotVarNames[i] == oxName) idOin = i;
      if (plotVarNames[i] == prodName) idPin = i;
      if (plotVarNames[i] == PPRName) idPPRin = i;
    }
    if (idFin<0 || idOin<0 || idFPRin < 0 || idPPRin < 0 || idPin < 0)
      Abort("Cannot find required data in pltfile");
    const int nCompIn  = 5;
    int nCompOut = 7;
    
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    //in
    const int idFlocal = 0; // fuel in here
    const int idOlocal = 1; // O2 in here
    const int idPlocal = 2; //product in here
    const int idFPRlocal = 3; //FCR in here
    const int idPPRlocal = 4; //Temp in here
    //out
    const int idZlocal = 0; // Z out here
    const int idCFlocal = 1; // CF out here
    const int idCPlocal = 2; // CP out here
    const int idCFsourcelocal = 3; //CF source here
    const int idCPsourcelocal = 4; //CP source here
    const int idCFCPlocal = 5; //(1-CF)*CP
    const int idCPCFlocal = 6; //(1-CP)*CF
    
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
    }
    inNames[idFlocal] = spName;
    inNames[idOlocal] =  oxName;
    inNames[idFPRlocal] = FPRName;
    inNames[idPPRlocal] = PPRName;
    inNames[idPlocal] = prodName;
    outNames[idZlocal] = "mixture_fraction";
    outNames[idCFlocal] = "progress_variable_F";
    outNames[idCPlocal] = "progress_variable_P";
    outNames[idCFsourcelocal] = "progress_variable_F_source";
    outNames[idCPsourcelocal] = "progress_variable_P_source";
    outNames[idCFCPlocal] = "CP-CPCF";
    outNames[idCPCFlocal] = "CF-CPCF";
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<MultiFab*> indata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    const int nGrow = 1;
    	  
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
	  Real CPfactor = Z_st+Y_O_air*(1-Z_st);
	  ZC(i,j,k,idZlocal) = (s*inbox(i,j,k,idFlocal) - inbox(i,j,k,idOlocal) + Y_O_air)/(s+Y_O_air);
	  if (ZC(i,j,k,idZlocal) <= Z_st){
	      ZC(i,j,k,idCFlocal) = 1-inbox(i,j,k,idFlocal)/ZC(i,j,k,idZlocal);
	      ZC(i,j,k,idCFsourcelocal) = -inbox(i,j,k,idFPRlocal)/ZC(i,j,k,idZlocal);
	      CPfactor *= ZC(i,j,k,idZlocal)/Z_st;
	  } else {
	      ZC(i,j,k,idCFlocal) = ((1-Z_st)/Z_st)*(ZC(i,j,k,idZlocal)-inbox(i,j,k,idFlocal))/(1-ZC(i,j,k,idZlocal));
	      ZC(i,j,k,idCFsourcelocal) = -((1-Z_st)/Z_st)*inbox(i,j,k,idFPRlocal)/(1-ZC(i,j,k,idZlocal));
	      CPfactor *= (1-ZC(i,j,k,idZlocal))/(1-Z_st);
	  }
	  ZC(i,j,k,idCPlocal) = inbox(i,j,k,idPlocal)/CPfactor;
	  ZC(i,j,k,idCFsourcelocal) = inbox(i,j,k,idPPRlocal)/CPfactor;
	  if (clipProgress){ 
	    ZC(i,j,k,idCFlocal) = std::max(0.0,ZC(i,j,k,idCFlocal));
	    ZC(i,j,k,idCFlocal) = std::min(1.0,ZC(i,j,k,idCFlocal));
	    ZC(i,j,k,idCPlocal) = std::max(0.0,ZC(i,j,k,idCPlocal));
	    ZC(i,j,k,idCPlocal) = std::min(1.0,ZC(i,j,k,idCPlocal));
	  }
	  ZC(i,j,k,idCFCPlocal) = (1-ZC(i,j,k,idCFlocal))*ZC(i,j,k,idCPlocal);
	  ZC(i,j,k,idCPCFlocal) = (1-ZC(i,j,k,idCPlocal))*ZC(i,j,k,idCFlocal);
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
