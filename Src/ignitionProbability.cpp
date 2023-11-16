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

//using namespace analysis_util;


Real closestIDT(Real H2, Real O2, Real T, Vector<Real> H2table, Vector<Real> O2table, Vector<Real> temptable, Vector<Real> IDTtable, int tableLength) {
  Vector<Real> stateDist(tableLength);
  for (int i = 0; i < tableLength; i++) {
    stateDist[i] = std::abs(H2-H2table[i])/(std::abs(H2)+std::abs(H2table[i])) + std::abs(O2-O2table[i])/(std::abs(O2)+std::abs(O2table[i])) + std::abs(T-temptable[i])/(std::abs(T)+std::abs(temptable[i]));
  }
  int index = std::min_element(stateDist.begin(),stateDist.end())-stateDist.begin();
  return IDTtable[index];
}

Real closestTTI(Real H2, Real O2, Real T, Real HO2, Vector<Real> H2table, Vector<Real> O2table, Vector<Real> temptable, Vector<Real> HO2table, Vector<Real> TTItable, int tableLength) {
  Vector<Real> stateDist(tableLength);
  for (int i = 0; i < tableLength; i++) {
    stateDist[i] = std::abs(H2-H2table[i])/(std::abs(H2)+std::abs(H2table[i])) + std::abs(O2-O2table[i])/(std::abs(O2)+std::abs(O2table[i])) + std::abs(T-temptable[i])/(std::abs(T)+std::abs(temptable[i]))+ std::abs(HO2-HO2table[i])/(std::abs(HO2)+std::abs(HO2table[i]));
  }
  int index = std::min_element(stateDist.begin(),stateDist.end())-stateDist.begin();
  return TTItable[index];
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
    Real tol = 1e-7; pp.query("tol",tol);
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
    Vector<std::string> varNames={"Y(H2)","Y(O2)","temp","Y(HO2)"};
    
    const int nCompIn  = 4;
    const int nCompOut = 4;
    Vector<std::string> outNames(nCompOut);
    Vector<int> destFillComps(nCompIn);
    
    for (int i = 0; i < nCompIn; i++) {
      destFillComps[i] = i;
    }
    outNames[0] = "IDT";
    outNames[1] = "p(I|IDT)";
    outNames[2] = "TTI";
    outNames[3] = "p(I|TTI)";
    //read in cantera data
    FILE *fp;
    Real input;

    fp = fopen("IDT.dat","r");
    Vector<Real> IDT;
    while(fscanf(fp,"%lf ", &input) == 1) {
      IDT.push_back(input);
    }
    fp = fopen("H2IDT.dat","r");
    Vector<Real> H2IDT;
    while(fscanf(fp,"%lf ", &input) == 1) {
      H2IDT.push_back(input);
    }
    fp = fopen("O2IDT.dat","r");
    Vector<Real> O2IDT;
    while(fscanf(fp,"%lf ", &input) == 1) {
      O2IDT.push_back(input);
    }
    fp = fopen("TempIDT.dat","r");
    Vector<Real> tempIDT;
    while(fscanf(fp,"%lf ", &input) == 1) {
      tempIDT.push_back(input);
    }
    //TTI
    fp = fopen("TTI.dat","r");
    Vector<Real> TTI;
    while(fscanf(fp,"%lf ", &input) == 1) {
      TTI.push_back(input);
    }
    fp = fopen("H2TTI.dat","r");
    Vector<Real> H2TTI;
    while(fscanf(fp,"%lf ", &input) == 1) {
      H2TTI.push_back(input);
    }
    fp = fopen("O2TTI.dat","r");
    Vector<Real> O2TTI;
    while(fscanf(fp,"%lf ", &input) == 1) {
      O2TTI.push_back(input);
    }
    fp = fopen("TempTTI.dat","r");
    Vector<Real> tempTTI;
    while(fscanf(fp,"%lf ", &input) == 1) {
      tempTTI.push_back(input);
    }
    fp = fopen("HO2TTI.dat","r");
    Vector<Real> HO2TTI;
    while(fscanf(fp,"%lf ", &input) == 1) {
      HO2TTI.push_back(input);
    }
    

    
    int tableLengthIDT = IDT.size();
    int tableLengthTTI = TTI.size();
    //for (int i = 0; i < tableLength; i++) {
    //  std::cout << IDT[i] << " " << H2IDT[i] << " " << O2IDT[i] << " " << tempIDT[i] << std::endl; 
    //}
   
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
      outdata[lev]->setVal(-1);
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
	  if  (inbox(i,j,k,0) < tol || inbox(i,j,k,1) < tol || inbox(i,j,k,2) < 800) {
	    outbox(i,j,k,0) = 1;
	    outbox(i,j,k,1) = 0;
	    outbox(i,j,k,2) = 1;
	    outbox(i,j,k,3) = 0;
	  } else {
	    outbox(i,j,k,0) = closestIDT(inbox(i,j,k,0),inbox(i,j,k,1),inbox(i,j,k,2),H2IDT,O2IDT,tempIDT,IDT,tableLengthIDT);
	    outbox(i,j,k,1) = std::exp(-50*outbox(i,j,k,0));
	    outbox(i,j,k,2) = closestTTI(inbox(i,j,k,0),inbox(i,j,k,1),inbox(i,j,k,2),inbox(i,j,k,3),H2TTI,O2TTI,tempTTI,HO2TTI,TTI,tableLengthTTI);
	    outbox(i,j,k,3) = std::exp(-50*outbox(i,j,k,2));
	  }
	});
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_IP");
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
