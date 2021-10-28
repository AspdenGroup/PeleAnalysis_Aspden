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
/*
#include <mechanism.h>
#include <chemistry_file.H>
#include <util.H>
#include <util_F.H>
#include <Transport_F.H>
#include <Fuego_EOS.H>
*/
using namespace amrex;
//using namespace analysis_util;
extern "C" {
    void pushvtog(const int* lo,  const int* hi,
                  const int* dlo, const int* dhi,
                  Real* U, const int* Ulo, const int* Uhi,
                  const int* nc);
    void pushvtog2d(const int* lo,  const int* hi,
                  const int* dlo, const int* dhi,
                  Real* U, const int* Ulo, const int* Uhi,
                  const int* nc);
}


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
      AmrData::SetVerbose(false);

    std::string plotFileName; pp.get("infile",plotFileName);
    std::string fuelName = "H2"; pp.query("fuelName",fuelName);
    Real s = 8.0; pp.query("stoichRatio",s);
    Real Y_O_air = 0.233; pp.query("YO2Air",Y_O_air);
    const Real Z_st = Y_O_air/(s+Y_O_air);
    Real D = 1.0; pp.query("diffCoeff",D);
    Real Le = 0.3; pp.query("Le",Le);
    int plotTerms = 0; pp.query("plotTerms",plotTerms);
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
    //Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string FCRName = fuelName+"_ConsumptionRate";
    const std::string spName= "Y("+fuelName+")";
    const std::string oxName= "Y(O2)";
    const std::string rhoName = "density";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == FCRName) idFCRin = i;
      if (plotVarNames[i] == spName) idFin = i;
      if (plotVarNames[i] == oxName) idOin = i;
      if (plotVarNames[i] == rhoName) idRhoin = i;
    }
    if (idFin<0 || idOin<0 || idFCRin < 0 || idRhoin < 0)
      Print() << "Cannot find required data in pltfile" << std::endl;
    const int nCompIn  = 4;
    int nCompOut = 3;
    if (plotTerms == 1) {
      nCompOut += 16;
    }
    
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    //in
    const int idFlocal = 0; // fuel in here
    const int idOlocal = 1; // O2 in here
    const int idRholocal = 2; // density in here
    const int idFCRlocal = 3; //FCR in here
    
    //out
    const int idZlocal = 0; // Z out here
    const int idClocal = 1; // C out here
    int idCdiff1local = -1, idZdiff1local = -1, idCSdiss1local = -1;
    int idZdiff2local = -1, idSdisslocal = -1, idCSdiss2local = -1, idCdiff2local = -1;
    int idZdiff2Clocal = -1, idSdissClocal = -1, idCSdiss2Clocal = -1, idCdiff2Clocal = -1;
    int idLewislocal = -1, idLewisClocal = -1;
    int idCsourcelocal = 2;
    if (plotTerms == 1) {
      idCdiff1local = 2; // C diffusion 1 term (standard)
      idZdiff1local = 3; // Z diffusion 1 term (standard)
      idCSdiss1local = 4; //Cross scalar dissipation 1 term
      //Z equation terms
      idZdiff2local = 5; //Z diffusion 2 term
      idSdisslocal = 6; //scalar dissipation term
      idCSdiss2local = 7; //cross scalar dissipation 2 term
      idCdiff2local = 8; //C diffusion 2 term
      //C equation terms
      idZdiff2Clocal = 9; //Z diffusion 2 term
      idSdissClocal = 10; //scalar dissipation term
      idCSdiss2Clocal = 11; //cross scalar dissipation 2 term
      idCdiff2Clocal = 12; //C diffusion 2 term
      
      idLewislocal = 13; //Lewis number term Z
      idLewisClocal = 14; //Lewis number term C
      idCsourcelocal += 13; //omega_C term
    }

    destFillComps[idFlocal] = idFlocal;
    destFillComps[idOlocal] = idOlocal;
    destFillComps[idRholocal] = idRholocal;
    destFillComps[idFCRlocal] = idFCRlocal;
    
    inNames[idFlocal] = spName;
    inNames[idOlocal] =  oxName;
    inNames[idRholocal] = rhoName;
    inNames[idFCRlocal] = FCRName;
    outNames[idZlocal] = "mixture_fraction";
    outNames[idClocal] = "progress_variable";
    if (plotTerms == 1) {
      outNames[idCdiff1local] = "standard_progress_variable_diffusion";
      outNames[idZdiff1local] = "standard_mixture_fraction_diffusion";
      outNames[idCSdiss1local] = "cross_scalar_dissipation_term_1";
      //Z equation terms
      outNames[idZdiff2local] = "mixture_fraction_diffusion_term_Z";
      outNames[idSdisslocal] = "scalar_dissipation_term_Z";
      outNames[idCSdiss2local] = "cross_scalar_dissipation_term_2_Z";
      outNames[idCdiff2local] = "progress_variable_diffusion_term_Z";
      //C equation terms
      outNames[idZdiff2Clocal] = "mixture_fraction_diffusion_term_C";
      outNames[idSdissClocal] = "scalar_dissipation_term_C";
      outNames[idCSdiss2Clocal] = "cross_scalar_dissipation_term_2_C";
      outNames[idCdiff2Clocal] = "progress_variable_diffusion_term_C";
      
      outNames[idLewislocal] = "Lewis_number_term_Z";
      outNames[idLewisClocal] = "Lewis_number_term_C";
    }
    outNames[idCsourcelocal] = "progress_variable_source";
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<MultiFab*> baseComps(Nlev);
    Vector<MultiFab*> indata(Nlev);
    Vector<Geometry*> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    const int nGrow = 1;
    int b[3] = {1, 1, 1};
    const int nBaseComps = 2;
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const Vector<Real>& delta = amrData.DxLevel()[lev];
      const DistributionMapping dm(ba);
      indata[lev] = new MultiFab(ba,dm,nCompIn,nGrow);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      //baseComps[lev] = new MultiFab(ba,dm,nBaseComps,nGrow);
      //MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(*indata[lev],lev,inNames,destFillComps); //Problem
      geoms[lev] = new Geometry(amrData.ProbDomain()[lev],&rb,0,&(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
        Array4<Real> const& ZC = (*outdata[lev]).array(mfi);
        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  ZC(i,j,k,idZlocal) = (s*inbox(i,j,k,idFlocal) - inbox(i,j,k,idOlocal) + Y_O_air)/(s+Y_O_air);
	  if (ZC(i,j,k,idZlocal) <= Z_st) {
	    ZC(i,j,k,idClocal) = 1-inbox(i,j,k,idFlocal)/ZC(i,j,k,idZlocal);
	    ZC(i,j,k,idCsourcelocal) = inbox(i,j,k,idFCRlocal)/ZC(i,j,k,idZlocal);
	  } else {
	    ZC(i,j,k,idClocal) = ((1-Z_st)/Z_st)*(ZC(i,j,k,idZlocal)-inbox(i,j,k,idFlocal))/(1-ZC(i,j,k,idZlocal));
	    ZC(i,j,k,idCsourcelocal) = ((1-Z_st)/Z_st)*inbox(i,j,k,idFCRlocal)/(1-ZC(i,j,k,idZlocal));
	  }
        });
	
      }

      if (plotTerms == 1) {
      
	//EXTRA STUFF TO FIX CF INTERFACE, NEEDS F90
	const Box& dbox = amrData.ProbDomain()[lev];
	for (MFIter mfi(*outdata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	  FArrayBox& fab = (*outdata[lev])[mfi];
	  const Box& box = mfi.validbox();
#if (BL_SPACEDIM == 2)
	  pushvtog2d(BL_TO_FORTRAN_BOX(box),BL_TO_FORTRAN_BOX(dbox),BL_TO_FORTRAN_ANYD(fab),&nBaseComps);
#else
	  pushvtog(BL_TO_FORTRAN_BOX(box),BL_TO_FORTRAN_BOX(dbox),BL_TO_FORTRAN_ANYD(fab),&nBaseComps);
#endif
	
	}
	outdata[lev]->FillBoundary(0,nBaseComps,geoms[lev]->periodicity());      
	/////
	//derivatives etc.
	
	for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
	  {
	    
	    const Box& bx = mfi.tilebox();
	    Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
	    Array4<Real> const& ZC = (*outdata[lev]).array(mfi);
	    Real tmpx,tmpy;
	    Real tmpz = 0;
	    Real Cdiff,Zdiff,CSdiss,Sdiss;
	    AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
				    {
				      //progress variable diffusion
				      tmpx = (inbox(i+1,j,k,idRholocal)-inbox(i+1,j,k,idRholocal))*(ZC(i+1,j,k,idClocal) - ZC(i-1,j,k,idClocal))/(4*delta[0]*delta[0]);
				      tmpx += inbox(i,j,k,idRholocal)*(ZC(i+1,j,k,idClocal) + ZC(i-1,j,k,idClocal)-2*ZC(i,j,k,idClocal))/(delta[0]*delta[0]);
				      
				      tmpy = (inbox(i,j+1,k,idRholocal)-inbox(i,j+1,k,idRholocal))*(ZC(i,j+1,k,idClocal) - ZC(i,j+1,k,idClocal))/(4*delta[1]*delta[1]);
				      tmpy += inbox(i,j,k,idRholocal)*(ZC(i,j+1,k,idClocal) + ZC(i,j-1,k,idClocal)-2*ZC(i,j,k,idClocal))/(delta[1]*delta[1]);
#if (BL_SPACEDIM == 3)
				      tmpz = (inbox(i,j,k+1,idRholocal)-inbox(i,j,k-1,idRholocal))*(ZC(i,j,k+1,idClocal) - ZC(i,j,k-1,idClocal))/(4*delta[2]*delta[2]);
				      tmpz += inbox(i,j,k,idRholocal)*(ZC(i,j,k+1,idClocal) + ZC(i,j,k-1,idClocal)-2*ZC(i,j,k,idClocal))/(delta[2]*delta[2]);
#endif
				      Cdiff = D*(tmpx+tmpy+tmpz);
				      //mixture fraction diffusion
				      tmpx = (inbox(i+1,j,k,idRholocal)-inbox(i+1,j,k,idRholocal))*(ZC(i+1,j,k,idZlocal) - ZC(i-1,j,k,idZlocal))/(4*delta[0]*delta[0]);
				      tmpx += inbox(i,j,k,idRholocal)*(ZC(i+1,j,k,idZlocal) + ZC(i-1,j,k,idZlocal)-2*ZC(i,j,k,idZlocal))/(delta[0]*delta[0]);
				      
				      tmpy = (inbox(i,j+1,k,idRholocal)-inbox(i,j+1,k,idRholocal))*(ZC(i,j+1,k,idClocal) - ZC(i,j+1,k,idClocal))/(4*delta[1]*delta[1]);
				      tmpy += inbox(i,j,k,idRholocal)*(ZC(i,j+1,k,idZlocal) + ZC(i,j-1,k,idZlocal)-2*ZC(i,j,k,idZlocal))/(delta[1]*delta[1]);
#if (BL_SPACEDIM == 3)
				      tmpz = (inbox(i,j,k+1,idRholocal)-inbox(i,j,k-1,idRholocal))*(ZC(i,j,k+1,idClocal) - ZC(i,j,k-1,idClocal))/(4*delta[2]*delta[2]);
				      tmpz += inbox(i,j,k,idRholocal)*(ZC(i,j,k+1,idZlocal) + ZC(i,j,k-1,idZlocal)-2*ZC(i,j,k,idZlocal))/(delta[2]*delta[2]);
#endif
				      Zdiff = D*(tmpx+tmpy+tmpz);
				      //scalar dissipation
				      tmpx = (ZC(i+1,j,k,idZlocal)-ZC(i-1,j,k,idZlocal))*(ZC(i+1,j,k,idZlocal)-ZC(i-1,j,k,idZlocal))/(4*delta[0]*delta[0]);
				      
				      tmpy = (ZC(i,j+1,k,idZlocal)-ZC(i,j-1,k,idZlocal))*(ZC(i,j+1,k,idZlocal)-ZC(i,j-1,k,idZlocal))/(4*delta[1]*delta[1]);
#if (BL_SPACEDIM == 3)	    
				      tmpz = (ZC(i,j,k+1,idZlocal)-ZC(i,j,k-1,idZlocal))*(ZC(i,j,k+1,idZlocal)-ZC(i,j,k-1,idZlocal))/(4*delta[2]*delta[2]);
#endif
				      Sdiss = 2*D*(tmpx+tmpy+tmpz);
				      //cross scalar dissipation
				      tmpx = (ZC(i+1,j,k,idClocal)-ZC(i-1,j,k,idClocal))*(ZC(i+1,j,k,idZlocal)-ZC(i-1,j,k,idZlocal))/(4*delta[0]*delta[0]);
				      
				      tmpy = (ZC(i,j+1,k,idClocal)-ZC(i,j-1,k,idClocal))*(ZC(i,j+1,k,idZlocal)-ZC(i,j-1,k,idZlocal))/(4*delta[1]*delta[1]);
#if (BL_SPACEDIM == 3)	    
				      tmpz = (ZC(i,j,k+1,idClocal)-ZC(i,j,k-1,idClocal))*(ZC(i,j,k+1,idZlocal)-ZC(i,j,k-1,idZlocal))/(4*delta[2]*delta[2]);
#endif
				      CSdiss = 2*D*(tmpx+tmpy+tmpz);
				      
				      ZC(i,j,k,idCdiff1local) = Cdiff;
				      ZC(i,j,k,idZdiff1local) = Zdiff;
				      ZC(i,j,k,idCSdiss1local) = inbox(i,j,k,idRhoin)*CSdiss/ZC(i,j,k,idZlocal);
				      //Z equation terms
				      ZC(i,j,k,idZdiff2local) = (1/Le - 1)*(1-ZC(i,j,k,idZlocal))*(1-ZC(i,j,k,idClocal))*Zdiff;
				      ZC(i,j,k,idSdisslocal) = -(1/Le -1)*(1-ZC(i,j,k,idClocal))*inbox(i,j,k,idRholocal)*Sdiss/2.0;
				      ZC(i,j,k,idCSdiss2local) = -(1/Le - 1)*(2-3*ZC(i,j,k,idZlocal))*inbox(i,j,k,idRholocal)*CSdiss/2.0;
				      ZC(i,j,k,idCdiff2local) = -(1/Le - 1)*ZC(i,j,k,idZlocal)*(1-ZC(i,j,k,idZlocal))*Cdiff;
				      //C equation terms
				      ZC(i,j,k,idZdiff2Clocal) = ZC(i,j,k,idZdiff1local)*(1-ZC(i,j,k,idClocal))/ZC(i,j,k,idZlocal);
				      ZC(i,j,k,idSdissClocal) = ZC(i,j,k,idSdisslocal)*(1-ZC(i,j,k,idClocal))/ZC(i,j,k,idZlocal);
				      ZC(i,j,k,idCSdiss2Clocal) = ZC(i,j,k,idCSdiss2local)*(1-ZC(i,j,k,idClocal))/ZC(i,j,k,idZlocal);
				      ZC(i,j,k,idCdiff2Clocal) = ZC(i,j,k,idCdiff2local)*(1-ZC(i,j,k,idClocal))/ZC(i,j,k,idZlocal);
				      ZC(i,j,k,idLewislocal) = ZC(i,j,k,idZdiff2local)+ZC(i,j,k,idSdisslocal)+ZC(i,j,k,idCSdiss2local)+ZC(i,j,k,idCdiff2local);
				      ZC(i,j,k,idLewisClocal) = ZC(i,j,k,idZdiff2Clocal)+ZC(i,j,k,idSdissClocal)+ZC(i,j,k,idCSdiss2Clocal)+ZC(i,j,k,idCdiff2Clocal);
				      
				    });
	    
	  }
      }
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_ZC");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
