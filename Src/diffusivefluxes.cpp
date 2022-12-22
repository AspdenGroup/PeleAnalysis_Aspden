#include <string>
#include <iostream>

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
  std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> gradVar=<name>\n";
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
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    Vector<std::string> species;
    std::string infile        = "";  
    int finestLevel           = 1000;
    //int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.get("infile",infile);
    pp.query("finestLevel",finestLevel);
    int nSpecies = pp.countval("species");
    species.resize(nSpecies);
    pp.getarr("species",species);
    Vector<std::string> massFractions(nSpecies);
    Vector<std::string> diffusionCoeffs(nSpecies);
    Vector<std::string> reactionTerms(nSpecies);
    for (int n = 0; n<nSpecies; n++) {
      massFractions[n] = "Y("+species[n]+")";
      diffusionCoeffs[n] = "rhoD("+species[n]+")";
      reactionTerms[n] = "I_R("+species[n]+")";
    }
    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    // Gradient variable
    Vector<int> idY(nSpecies,-1);
    Vector<int> idD(nSpecies,-1);
    Vector<int> idR(nSpecies,-1);
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      for (int j=0; j < nSpecies; j++ ) {
	if (plotVarNames[i] == massFractions[j]) {
	  idY[j] = i;
	  break;
	}
	if (plotVarNames[i] == diffusionCoeffs[j]) {
	  idD[j] = i;
	  break;
	}
	if (plotVarNames[i] == reactionTerms[j]) {
	  idR[j] = i;
	  break;
	}
      }
    }
    for (int i=0; i<nSpecies; ++i) {
      if (idY[i]<0 || idD[i]<0 || idR[i] < 0) {
	Abort("Cannot find required data in pltfile");
      }
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    int nCompIn = 3*nSpecies;

    Vector<int> destFillComps(3*nSpecies);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }
    Vector<std::string> inVarNames(3*nSpecies);
    for (int i=0; i<nSpecies; i++) {
      inVarNames[i] = massFractions[i];
      inVarNames[i+nSpecies] = diffusionCoeffs[i];
      inVarNames[i+2*nSpecies] = reactionTerms[i];
    }
    
    const int idGr = 3*nSpecies;
    const int nCompGrads = idGr+AMREX_SPACEDIM*nSpecies; //Y's, D's, gradY's, all the work will be done here then copied into final MF
    const int nCompOut = 2*nSpecies; // Diffusive fluxes and Damkohler numbers 

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    BCRec gradVarBC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradVarBC.setLo(idim,BCType::foextrap);
        gradVarBC.setHi(idim,BCType::foextrap);
        if ( is_per[idim] ) {
            gradVarBC.setLo(idim, BCType::int_dir);
            gradVarBC.setHi(idim, BCType::int_dir);
        }
    }

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab> workingMF(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    Vector<MultiFab> outdata(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev=0; lev<Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      workingMF[lev].define(grids[lev], dmap[lev], nCompGrads, nGrow);
      outdata[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      
      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(workingMF[lev], lev, inVarNames, destFillComps);

      workingMF[lev].FillBoundary(0,1,geoms[lev].periodicity());
    }

    // Get face-centered gradients from MLMG
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);
    poisson.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson.setDomainBC(lo_bc, hi_bc);

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0;                   // No need for ghost face on gradient

    for (int n =0; n<nSpecies; n++) {
      Print() << "Calculating grad"+massFractions[n] << std::endl;
      Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
      Vector<std::unique_ptr<MultiFab>> phi;
      Vector<MultiFab> laps;
      for (int lev = 0; lev < Nlev; ++lev) {
	for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
	  const auto& ba = grids[lev];
	  grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
				 dmap[lev], 1, nGrowGrad);
	}    
	phi.push_back(std::make_unique<MultiFab> (workingMF[lev],amrex::make_alias,n,1));
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      Print() << "Averaging to cell centers and multiplying by "+diffusionCoeffs[n]<<std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
	MultiFab gradAlias(workingMF[lev], amrex::make_alias, idGr+n*AMREX_SPACEDIM, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
	for (int d = 0; d<AMREX_SPACEDIM; d++) {
	  MultiFab::Multiply(gradAlias,workingMF[lev],nSpecies+n,d,1,0);
	}
      }
    } //we have rhoDgradY, differentiate again
    
    for (int n =0; n<AMREX_SPACEDIM*nSpecies; n++) {
      int d = n%AMREX_SPACEDIM; //keep track of which derivative
      Print() << "Differentiating " +diffusionCoeffs[n/AMREX_SPACEDIM] + "grad"+massFractions[n/AMREX_SPACEDIM] << std::endl;
      Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
      Vector<std::unique_ptr<MultiFab>> phi;
      Vector<MultiFab> laps;
      for (int lev = 0; lev < Nlev; ++lev) {
	for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
	  const auto& ba = grids[lev];
	  grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
				 dmap[lev], 1, nGrowGrad);
	}    
	phi.push_back(std::make_unique<MultiFab> (workingMF[lev],amrex::make_alias,idGr+n,1));
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      Print() << "Averaging to cell centers and only grabbing derivative in direction "+std::to_string(d) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
	Print() << "Making CC MF\n";
	MultiFab gradCC(grids[lev],dmap[lev],AMREX_SPACEDIM,nGrow); //need to store all gradients for averaging, temporarily
	Print() << "Averaging to CC\n";
	average_face_to_cellcenter(gradCC, 0, GetArrOfConstPtrs(grad[lev]));
	Print() << "Copying back to workingMF\n";
	MultiFab::Copy(workingMF[lev],gradCC,d,idGr+n,1,nGrow); //just take dim we want
      }
    }
    Print() << "Building output data...\n";
    //time to build MF with data to output
    for (int lev = 0; lev<Nlev; ++lev) {
      outdata[lev].setVal(0);
      for (int n = 0; n < nSpecies; n++) {
	for (int d = 0; d < AMREX_SPACEDIM; d++) {
	  MultiFab::Add(outdata[lev],workingMF[lev],idGr+n*AMREX_SPACEDIM+d,n,1,nGrow); //add x,y and z derivatives
	}
      }
      MultiFab::Copy(outdata[lev],workingMF[lev],2*nSpecies,nSpecies,nSpecies,nGrow); //copy reactions in
      MultiFab::Divide(outdata[lev],outdata[lev],0,nSpecies,nSpecies,nGrow); //divide reactions by diffusive flux
      outdata[lev].abs(nSpecies, nSpecies); //take absolute value
    }
    
    
    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    for (int i=0; i<nSpecies; ++i) {
      nnames[i] = "-grad("+diffusionCoeffs[i]+"grad"+massFractions[i]+")";
      nnames[i+nSpecies] = "Da("+species[i]+")";
    }
    std::string outfile(getFileRoot(infile) + "_DF"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), nnames,
                                   geoms, 0.0, isteps, refRatios);
  }
  amrex::Finalize();
  return 0;
}
