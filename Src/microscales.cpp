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
    
    int idCSin = -1;
    int idSijin = -1;
    int idRhoin = -1;
    int idMuin = -1;
    int idDFin = -1;
    int idDOin = -1;
    //Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string cSourceName = "progress_variable_source";
    const std::string sijName = "SijSij";
    const std::string rhoName = "density";
    const std::string muName = "mu";
    const std::string DFuelName = "rhoD("+fuelName+")";
    const std::string DOxName = "rhoD(O2)";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == cSourceName) idCSin = i;
      if (plotVarNames[i] == sijName) idSijin = i;
      if (plotVarNames[i] == rhoName) idRhoin = i;
      if (plotVarNames[i] == muName) idMuin = i;
      if (plotVarNames[i] == DFuelName) idDFin = i;
      if (plotVarNames[i] == DOxName) idDOin = i;
    }
    if (idCSin < 0 || idSijin < 0 || idRhoin < 0 || idMuin < 0 || idDFin < 0 || idDOin < 0)
      Abort("Cannot find required data in pltfile");
    const int nCompIn  = 6;
    int nCompOut = 10; //tau_c, tau_eta, Da, eps, eta, eta_H2, eta_O2
    
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    //in
    const int idCSlocal = 0; 
    const int idSijlocal = 1;
    const int idRholocal = 2;
    const int idMulocal = 3; 
    const int idDFlocal = 4; 
    const int idDOlocal = 5;
    //out
    const int idTClocal = 0; 
    const int idTElocal = 1; 
    const int idDalocal = 2; 
    const int idEpslocal = 3; 
    const int idEtalocal = 4;
    const int idEtaFlocal = 5;
    const int idEtaOlocal = 6; 
    const int idEtaDxlocal = 7;
    const int idEtaFDxlocal = 8;
    const int idEtaODxlocal = 9;
    for (int n = 0; n<nCompIn; n++) {
      destFillComps[n] = n;
    }
    inNames[idCSlocal] = cSourceName;
    inNames[idSijlocal] =  sijName;
    inNames[idRholocal] = rhoName;
    inNames[idMulocal] = muName;
    inNames[idDFlocal] = DFuelName;
    inNames[idDOlocal] = DOxName;
    outNames[idTClocal] = "tau_c";
    outNames[idTElocal] = "tau_eta";
    outNames[idDalocal] = "log(Da)";
    outNames[idEpslocal] = "dissipation_rate";
    outNames[idEtalocal] = "eta";
    outNames[idEtaFlocal] = "eta_"+fuelName;
    outNames[idEtaOlocal] = "eta_O2";
    outNames[idEtaDxlocal] = "eta_dx";
    outNames[idEtaFDxlocal] = "eta_"+fuelName+"_dx";
    outNames[idEtaODxlocal] = "eta_O2_dx";
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<MultiFab*> indata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    const int nGrow = 0;
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
	Array4<Real> const& inArr  = (*indata[lev]).array(mfi);
        Array4<Real> const& outArr = (*outdata[lev]).array(mfi);
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  Real RR = std::max(1e-12, inArr(i,j,k,idCSlocal));
	  outArr(i,j,k,idTClocal) = 1.0/RR;
	  outArr(i,j,k,idTElocal) = 1.0/std::sqrt(2.0*inArr(i,j,k,idSijlocal));
	  outArr(i,j,k,idDalocal) = std::log10(outArr(i,j,k,idTElocal)/outArr(i,j,k,idTClocal));
	  Real kin_visc = inArr(i,j,k,idMulocal)/inArr(i,j,k,idRholocal);
	  outArr(i,j,k,idEpslocal) = std::max(1e-12,2*kin_visc*inArr(i,j,k,idSijlocal));
	  outArr(i,j,k,idEtalocal) = std::pow(std::pow(kin_visc,3.0)/outArr(i,j,k,idEpslocal),0.25);
	  Real DF = inArr(i,j,k,idDFlocal)/inArr(i,j,k,idRholocal); 
	  outArr(i,j,k,idEtaFlocal) =  std::pow(std::pow(DF,3.0)/outArr(i,j,k,idEpslocal),0.25);
	  Real DO = inArr(i,j,k,idDOlocal)/inArr(i,j,k,idRholocal); 
	  outArr(i,j,k,idEtaOlocal) =  std::pow(std::pow(DO,3.0)/outArr(i,j,k,idEpslocal),0.25);
	  outArr(i,j,k,idEtaDxlocal) = outArr(i,j,k,idEtalocal)/delta[0];
	  outArr(i,j,k,idEtaFDxlocal) = outArr(i,j,k,idEtaFlocal)/delta[0];
	  outArr(i,j,k,idEtaODxlocal) = outArr(i,j,k,idEtaOlocal)/delta[0];
	});
	
      }
      
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_microscales");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
