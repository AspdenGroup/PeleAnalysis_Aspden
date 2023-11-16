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
    //int dualFuel = 0; pp.query("duelFuel");
    std::string fuelName = "H2"; pp.query("fuelName",fuelName);
    //if (dualFuel) {
    std::string fuelName2=""; pp.query("fuelName2",fuelName2);
    //}
    int fuelWeight = 1; pp.query("fuelWeight",fuelWeight);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    bool initNormalise=false;
    int initQuery = 1; pp.query("initNormalise",initQuery); 

    DataServices dataServicesInit("plt00000", fileType);
    if (!dataServicesInit.AmrDataOk()) {
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Cannot find initial condition - normalising using current plotfile" << std::endl;
    } else if (initQuery == 0){
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Manually set to not normalise" << std::endl;
    } else {
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Normalising using initial condition" << std::endl;
      initNormalise=true;
    }

    

    AmrData& initAmrData = dataServicesInit.AmrDataRef();
    

    //init_mech();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    int idTin = -1;
    //Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y("+fuelName+")";
    
    //if (dualFuel) {
    const std::string spName2 = "Y("+fuelName2+")";
      //}
    const std::string TName= "temp";
    //std::cout << spName << std::endl;
    /*for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName) idTin = i;
    }
    if (idYin<0 || idTin<0)
    Print() << "Cannot find required data in pltfile" << std::endl;*/
    const int dualFuel = !fuelName2.empty(); 
    int nCompIn = 2;
    if (dualFuel) {
      nCompIn += 1;
    }
    const int nCompOut = 2;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    const int idTlocal = 0; // T here
    const int idYlocal = 1; // Y here
    const int idY2local = 2;
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idYlocal] = idYlocal;
    if (dualFuel) {
      destFillComps[idY2local] = idY2local;
    }
    
    inNames[idTlocal] = TName;
    inNames[idYlocal] =  spName;
    outNames[idTlocal] = "prog_temp";
    if (dualFuel) {
      inNames[idY2local] = spName2;
      outNames[idYlocal] = "prog_blend";
    } else {
      outNames[idYlocal] = "prog_"+fuelName;
    }
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;
    int b[3] = {1, 1, 1};
    Real Y_fuel_u, Y_fuel_b, T_u, T_b, Y_blend_u, Y_blend_b;
    if(initNormalise) {
      if(initAmrData.MinMax(initAmrData.ProbDomain()[0],spName,0,Y_fuel_b,Y_fuel_u) && initAmrData.MinMax(initAmrData.ProbDomain()[0],TName,0,T_u,T_b)) {
	if (dualFuel && initAmrData.MinMax(initAmrData.ProbDomain()[0],spName2,0,Y_blend_b,Y_blend_u)) {
	  if (ParallelDescriptor::IOProcessor())
	  std::cout << "Found min/max (blend)" << std::endl;
	} else {
	  if (ParallelDescriptor::IOProcessor())
	    std::cout << "Found min/max" << std::endl;
	}
      } else {
	std::cout << "Could not find min/max" << std::endl;
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
      }
    } else {
       if(amrData.MinMax(amrData.ProbDomain()[0],spName,0,Y_fuel_b,Y_fuel_u) && amrData.MinMax(amrData.ProbDomain()[0],TName,0,T_u,T_b)){
	 if (dualFuel && amrData.MinMax(amrData.ProbDomain()[0],spName2,0,Y_blend_b,Y_blend_u)) {
	   if (ParallelDescriptor::IOProcessor())
	     std::cout << "Found min/max (blend)" << std::endl;
	} else {
	   if (ParallelDescriptor::IOProcessor())
	     std::cout << "Found min/max" << std::endl;
	 }
       } else {
	 std::cout << "Could not find min/max" << std::endl;
	 DataServices::Dispatch(DataServices::ExitRequest, NULL);
       }
    }

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //Problem
      Print() << "Data has been read for level " << lev << std::endl;
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& inbox  = indata.array(mfi);
        Array4<Real> const& outbox = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  outbox(i,j,k,idTlocal) = (inbox(i,j,k,idTlocal)-T_u)/(T_b-T_u);
	  if (dualFuel) {
	    outbox(i,j,k,idYlocal) = 1-(fuelWeight*inbox(i,j,k,idYlocal)+inbox(i,j,k,idY2local))/(fuelWeight*Y_fuel_u+Y_blend_u);
	  } else {
	    outbox(i,j,k,idYlocal) = 1-inbox(i,j,k,idYlocal)/Y_fuel_u;
	  }
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_prog");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    WritePlotFile(GetVecOfPtrs(outdata),amrData,outfile,verb,outNames);
  }
  Finalize();
  return 0;
}
