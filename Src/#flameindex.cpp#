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
    Vector<std::string> gradVars(2);
    std::string fuel="Y(H2)";
    std::string oxidizer="Y(O2)";
    std::string infile        = "";  
    int finestLevel           = 1000;
    int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.get("infile",infile);
    pp.query("finestLevel",finestLevel);
    pp.query("fuel",fuel);
    gradVars[0] = fuel;
    gradVars[1] = oxidizer;
    
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
    Vector<int> idVels(2);
    idVels[0] = -1;
    idVels[1] = -1;
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      for (int j=0; j < 2; j++ ) {
	if (plotVarNames[i] == gradVars[j]) {
	  idVels[j] = i;
	  break;
	}
      }
    }
    for (int i=0; i<2; ++i) {
      if (idVels[i]<0) {
	Print() << "Cannot find " << gradVars[i] << " data in pltfile \n";
      }
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    int nCompIn = 2;
    Vector<std::string> inVarNames(2);
    inVarNames = gradVars;

    Vector<int> destFillComps(2);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }
    
    const int idGr = nCompIn;
    const int nCompOut = idGr + AMREX_SPACEDIM*2 + 1; // All the gradients plus 

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
    Vector<MultiFab> state(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev=0; lev<Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);

      state[lev].FillBoundary(0,1,geoms[lev].periodicity());
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

    for (int n =0; n<2; n++) {
      Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
      Vector<std::unique_ptr<MultiFab>> phi;
      Vector<MultiFab> laps;
      for (int lev = 0; lev < Nlev; ++lev) {
	for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
	  const auto& ba = grids[lev];
	  grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
				 dmap[lev], 1, nGrowGrad);
	}    
	phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,n,1));
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      
      for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
	MultiFab gradAlias(state[lev], amrex::make_alias, idGr+n*AMREX_SPACEDIM, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      }
    }
    
    for (int lev = 0; lev<Nlev; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {    
           const Box& bx = mfi.tilebox();
           auto const& grad_a   = state[lev].const_array(mfi,idGr);
           auto const& FI  = state[lev].array(mfi,idGr+2*AMREX_SPACEDIM);
           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
	     FI(i,j,k) = AMREX_D_TERM(  grad_a(i,j,k,0) * grad_a(i,j,k,AMREX_SPACEDIM),
					     + grad_a(i,j,k,1) * grad_a(i,j,k,1+AMREX_SPACEDIM),
					      + grad_a(i,j,k,2) * grad_a(i,j,k,2+AMREX_SPACEDIM));
	     FI(i,j,k) /= std::sqrt(AMREX_D_TERM(grad_a(i,j,k,0) * grad_a(i,j,k,0), + 
						 grad_a(i,j,k,1) * grad_a(i,j,k,1), +
						 grad_a(i,j,k,2) * grad_a(i,j,k,2)));
	     FI(i,j,k) /= std::sqrt(AMREX_D_TERM(grad_a(i,j,k,AMREX_SPACEDIM) * grad_a(i,j,k,AMREX_SPACEDIM), + 
						 grad_a(i,j,k,1+AMREX_SPACEDIM) * grad_a(i,j,k,1+AMREX_SPACEDIM), +
						 grad_a(i,j,k,2+AMREX_SPACEDIM) * grad_a(i,j,k,2+AMREX_SPACEDIM)));
           });  
        }
    	
     }
    
    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }
    for (int i=0; i<2; i++) {
      nnames[idGr+i*AMREX_SPACEDIM+0] = gradVars[i] + "_gx";
      nnames[idGr+i*AMREX_SPACEDIM+1] = gradVars[i] + "_gy";
#if AMREX_SPACEDIM==3
      nnames[idGr+i*AMREX_SPACEDIM+2] = gradVars[i] + "_gz";
#endif
    }
    nnames[idGr+2*AMREX_SPACEDIM] = "FI";
    std::string outfile(getFileRoot(infile) + "_fi"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(state), nnames,
                                   geoms, 0.0, isteps, refRatios);
  }
  amrex::Finalize();
  return 0;
}
