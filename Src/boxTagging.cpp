#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

#define Y_IN_PLOTFILE
#undef Y_IN_PLOTFILE

static
void 
print_usage (int,
             char* argv[])
{
	std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> outfile=<> [options] \n\tOptions:\n";
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
    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    std::string infile;
    pp.get("infile",infile);

    std::string outfile;
    pp.get("outfile",outfile);

    int avgDown = 0; pp.query("avgDown",avgDown);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Set up for reading pltfile
    DataServices dataServices(infile, fileType);

    if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData& amrData = dataServices.AmrDataRef();

    int RRcomp = -1;
    int Tcomp = -1;
    Real Treac = 300;
    Real Tprod = 2000;
    Real RRthresholdC = 0;
    Real RRthresholdE = 0;
    pp.get("RRcomp",RRcomp);
    pp.get("RRthresholdC",RRthresholdC);
    pp.get("RRthresholdE",RRthresholdE);
    pp.get("Tcomp",Tcomp);
    pp.get("Treac",Treac);
    pp.get("Tprod",Tprod);
    
    int finestLevel = amrData.FinestLevel();
    Vector<Real> plo = amrData.ProbLo();
    Real productHeight = amrData.ProbHi()[2];
    pp.query("productHeight",productHeight);
    int Nlev = finestLevel + 1;
    pp.query("finestLevel",finestLevel);
    Nlev = std::max(0, std::min(Nlev, finestLevel+1));
    int max_grid_size = 8;
    //pp.query("max_grid_size",max_grid_size);
    Vector<MultiFab*> fileData(Nlev);
    int nGrowC = 1;
    int nGrowE = 1;
    pp.query("growSizeC",nGrowC);
    pp.query("growSizeE",nGrowE);
    //int nGrow = std::max(nGrowC,nGrowE);
    Vector<int> nGrowELev(Nlev);
    Vector<int> nGrowCLev(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      BoxArray newba = amrData.boxArray(lev);
      
      const DistributionMapping dm(newba);
      int nGrow = std::max(nGrowC,nGrowE)*std::pow(2,lev);
      nGrowELev[lev] = nGrowE*std::pow(2,lev);
      nGrowCLev[lev] = nGrowC*std::pow(2,lev);
      fileData[lev] = new MultiFab(newba,dm,3,nGrow);
      
      fileData[lev]->copy(amrData.GetGrids(lev,RRcomp),0,0,1);
      fileData[lev]->copy(amrData.GetGrids(lev,Tcomp),0,1,1);
      /*(if (lev == finestLevel) {
	smallBoxes[lev]->setVal(3,2,1);
      } else {
      */
      fileData[lev]->setVal(-1,2,1);
      Vector<Real> dx = amrData.DxLevel()[lev];
      for (MFIter mfi(*fileData[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	Array4<Real> const& varsIn = fileData[lev]->array(mfi);
	
	AMREX_PARALLEL_FOR_3D (bx,i,j,k, {
	    if(varsIn(i,j,k,2) < 0) {
	      Real z = plo[2]+(k+0.5)*dx[2];
	      if (varsIn(i,j,k,0) > RRthresholdC) { //Core
		for (int ii=-nGrowCLev[lev]; ii<=nGrowCLev[lev]; ii++) {
		  for (int jj=-nGrowCLev[lev]; jj<=nGrowCLev[lev]; jj++) {
		    for (int kk=-nGrowCLev[lev]; kk<=nGrowCLev[lev]; kk++) {
		      varsIn(i+ii,j+jj,k+kk,2) = 3;
		    }
		  }
		}
	      } else if (varsIn(i,j,k,1) >= Tprod || z > productHeight) { //Products
		varsIn(i,j,k,2) = 4;
	      } else if (varsIn(i,j,k,1) <= Treac) { //Jet
		varsIn(i,j,k,2) = 0;
	      } else if (varsIn(i,j,k,0) > RRthresholdE) { //Edge
		for (int ii=-nGrowELev[lev]; ii<=nGrowELev[lev]; ii++) {
		  for (int jj=-nGrowELev[lev]; jj<=nGrowELev[lev]; jj++) {
		    for (int kk=-nGrowELev[lev]; kk<=nGrowELev[lev]; kk++) {
		      varsIn(i+ii,j+jj,k+kk,2) = 2;
		    }
		  }
		}
	      } else { //Recirculation
		varsIn(i,j,k,2) = 1;
	      }
	    }
	  });	
      }
    }
    if (ParallelDescriptor::IOProcessor())
        std::cerr << "Set zones" << std::endl;
        
    //recast?
    Vector<std::string> names(1);
    names[0] = "zone";
    Vector<Geometry> geoms(Nlev);
    Vector<int> levelSteps(Nlev);
    Vector<IntVect> refRatios(Nlev-1);

    RealBox rb(&(amrData.ProbLo()[0]),
               &(amrData.ProbHi()[0]));
    Vector<int> is_per(BL_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,BL_SPACEDIM);
    int coord = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      if (lev < Nlev-1) {
        int r = amrData.RefRatio()[lev];
        refRatios[lev] = IntVect(D_DECL(r,r,r));
      }
    }
    if (avgDown) {
      Print() << "Averaging down to coarse and writing to "+outfile << std::endl;
      for (int lev = finestLevel; lev > 0; --lev) {
	MultiFab* fineMF = fileData[lev];
	MultiFab* coarseMF = fileData[lev-1];
	amrex::average_down(*fineMF,*coarseMF,0,1,refRatios[lev-1]);
      }
      amrex::WriteSingleLevelPlotfile(outfile,*fileData[0],names,geoms[0],0.0,0);
    } else {
      Print() << "Writing to " + outfile << std::endl;
      Vector<int> isteps(Nlev, 0);
      amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(fileData),names,geoms, 0.0, isteps, refRatios);
    }

    //WriteMultiLevelPlotfile(outfile,Nlev,dat,names,geoms,amrData.Time(),levelSteps,refRatio);
    }
    amrex::Finalize();
    return 0;
}
