#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

using namespace amrex;

Real interpolate(Real x, Real y, Real dx, Real dy, Vector<Vector<Real>> idttab) {
  Real xloc = x/dx;
  Real yloc = y/dy;
  int xloidx = (int)xloc;
  int yloidx = (int)yloc;
  int tabxsize = idttab.size();
  int tabysize = idttab[0].size();
  Real minidt = 1e-6;
  if (xloidx < 0 || xloidx+1 > tabxsize - 1 || yloidx < 0 || yloidx+1 > tabysize - 1) {
    return 0.0;
  } else if (idttab[xloidx][yloidx] < minidt || idttab[xloidx+1][yloidx] < minidt || idttab[xloidx][yloidx+1] < minidt ||idttab[xloidx+1][yloidx+1] < minidt) {
    return 0.0;
  } else {
    Real alphax = xloc-xloidx;
    Real alphay = yloc-yloidx;
    Real idtxlo = (1-alphax)*idttab[xloidx][yloidx] + alphax*idttab[xloidx+1][yloidx];
    Real idtxhi = (1-alphax)*idttab[xloidx][yloidx+1] + alphax*idttab[xloidx+1][yloidx+1];
    Real idt = (1-alphay)*idtxlo + alphay*idtxhi;
    return idt;
  }
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
    Vector<std::string> varNames={"Y(H2)","Y(H2O)"};
    
    const int nCompIn  = 2;
    const int nCompOut = 4;
    Vector<std::string> outNames(nCompOut);
    Vector<int> destFillComps(nCompIn);
    
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
    }
    outNames[0] = "x";
    outNames[1] = "y";
    outNames[2] = "IDT";
    outNames[3] = "1/IDT";
    //read in cantera data

    FILE *fp;
    Real input;
    
    fp = fopen("x.dat","r");
    Vector<Real> xtab;
    while(fscanf(fp,"%lf ", &input) == 1) {
      xtab.push_back(input);
    }
    fp = fopen("y.dat","r");
    Real dx = xtab[1]-xtab[0];
    Vector<Real> ytab;
    while(fscanf(fp,"%lf ", &input) == 1) {
      ytab.push_back(input);
    }
    Real dy = ytab[1]-ytab[0];
    int xsize = xtab.size();
    int ysize = ytab.size();
    
    std::ifstream fp2("IDT.dat");
    Vector<Real> tmp(ysize);
    Vector<Vector<Real>> idttab(xsize,tmp);
    if (!fp2) {
      Print() << "File couldn't be opened" << std::endl;
      return 1;
    }
    for (int row = 0; row < xsize; row++) {
      for (int column = 0; column < ysize; column++) {
	fp2 >> idttab[row][column];
	if (!fp2) {
	  Print() << "Error reading file" << std::endl;
	  return 1;
	}
      }
    }
    
    //Real Tair, Tfuel, Tb;
    //pp.get("Tfuel",Tfuel);
    //pp.get("Tair",Tair);
    //pp.get("Tb",Tb);
    Real H2Ob;//, H2u;
    pp.get("H2Ob",H2Ob);
    //pp.get("H2u",H2u);
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    const int nGrow = 0;
    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));
    int coord = 0;
    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,varNames,destFillComps); //Problem
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = indata.array(mfi);
        Array4<Real> const& outbox = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {	  
	  Real x = inbox(i,j,k,0);
	  Real y = inbox(i,j,k,1)/H2Ob;
	  outbox(i,j,k,0) = x;
	  outbox(i,j,k,1) = y;
	  Real idt = interpolate(x,y,dx,dy,idttab); //returns 0 if not autoignitable
	  if (idt > 1e-7) { 
	    outbox(i,j,k,2) = idt;
	    outbox(i,j,k,3) = 1/idt;
	  } else {  
	    outbox(i,j,k,2) = 0.01;
	    outbox(i,j,k,3) = 0;
	  }
	});
      }
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_IDT");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);

  }
  Finalize();
  return 0;
}
