#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <StreamPC.H>

using namespace amrex;

static
Vector<Vector<Real>>
GetSeedLocations (const StreamParticleContainer& spc, Vector<int>& faceData)
{
  Vector<Vector<Real>> locs;

  ParmParse pp;
  int nc=pp.countval("oneSeedPerCell");
  int ni=pp.countval("isofile");
  int ns=pp.countval("seedLoc");
  int nrL=pp.countval("seedRakeL");
  int nrR=pp.countval("seedRakeR");
  AMREX_ALWAYS_ASSERT((nc>0) ^ ((ni>0) ^ ((ns>0) ^ ((nrL>0) && nrR>0))));
  if (nc>0)
  {
    Abort();
    
    int finestLevel = spc.numLevels() - 1;
    std::vector< std::pair<int,Box> > isects;
    FArrayBox mask;
    for (int lev=0; lev<=finestLevel; ++lev)
    {
      const auto& geom = spc.Geom(lev);
      const auto& dx = geom.CellSize();
      const auto& plo = geom.ProbLo();

      BoxArray baf;
      if (lev < finestLevel) {
        baf = BoxArray(spc.ParticleBoxArray(lev+1)).coarsen(spc.GetParGDB()->refRatio(lev));
      }
      for (MFIter mfi = spc.MakeMFIter(lev); mfi.isValid(); ++mfi)
      {
        const Box& tile_box  = mfi.tilebox();
        if (AMREX_SPACEDIM<3 || tile_box.contains(IntVect(D_DECL(0,50,107)))) {

          mask.resize(tile_box,1);
          mask.setVal(1);
          if (lev < finestLevel) {
            isects = baf.intersections(tile_box);
            for (const auto& p : isects) {
              mask.setVal(0,p.second,0,1);
            }
          }

          for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
          {
            if (mask(iv,0) > 0)
            {
              locs.push_back({AMREX_D_DECL(plo[0] + (iv[0] + 0.5)*dx[0],
                                           plo[1] + (iv[1] + 0.5)*dx[1],
                                           plo[2] + (iv[2] + 0.5)*dx[2])});
            }
          }
        }
      }
    }
  }
  else if (ni>0)
  {
    // Read in isosurface
    AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM==3);
    std::string isofile; pp.get("isofile",isofile);
    Print() << "Reading isofile... " << isofile << std::endl;

    std::ifstream ifs;
    ifs.open(isofile.c_str());
    // AJA added dummy line read; sometimes time, sometimes `decimated'
    std::string topline;
    std::getline(ifs,topline);
    //Print() << "topline = " << topline << std::endl;
    std::string line;
    std::getline(ifs,line);
    auto surfNames = Tokenize(line,std::string(", "));
    int nCompSeedNodes = surfNames.size();
    int nElts, nodesPerElt;
    ifs >> nElts;
    ifs >> nodesPerElt;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    int nSeedNodes = tnodes.box().numPts();

    Real* ndat = tnodes.dataPtr();
    //Print() << "First few locs:" << std::endl;
    for (int i=0; i<nSeedNodes; ++i)
    {
      int o=i*nCompSeedNodes;
      locs.push_back({AMREX_D_DECL(ndat[o+0], ndat[o+1], ndat[o+2])});
      /*
      if (i<13) {
	Print() << locs[i][0] << " " << locs[i][1] << " " << locs[i][2] << std::endl;
      }
      */
    }
    tnodes.clear();

    faceData.resize(nElts*nodesPerElt);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ifs.close();
  }
  else if (pp.countval("seedLoc")>0)
  {
    Abort();
    Vector<Real> loc(AMREX_SPACEDIM);
    pp.getarr("seedLoc",loc,0,AMREX_SPACEDIM);
    locs.push_back({AMREX_D_DECL(loc[0], loc[1], loc[2])});
  }
  else
  {
    Abort();
    int seedRakeNum;
    pp.get("seedRakeNum",seedRakeNum);
    AMREX_ALWAYS_ASSERT(seedRakeNum >= 2);
    Vector<Real> locL(AMREX_SPACEDIM), locR(AMREX_SPACEDIM);
    pp.getarr("seedRakeL",locL,0,AMREX_SPACEDIM);
    pp.getarr("seedRakeR",locR,0,AMREX_SPACEDIM);

    for (int i=0; i<seedRakeNum; ++i) {
      locs.push_back({AMREX_D_DECL(locL[0] + (i/double(seedRakeNum-1))*(locR[0] - locL[0]),
                                   locL[1] + (i/double(seedRakeNum-1))*(locR[1] - locL[1]),
                                   locL[2] + (i/double(seedRakeNum-1))*(locR[2] - locL[2]))});
    }
  }  
  return locs;
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    // get base for output files
    std::string outfile = infile;
    pp.query("outfile",outfile);
    //
    int writeParticles(0);
    pp.query("writeParticles",writeParticles);
    std::string partfile = outfile+"_particles";
    //
    int writeStreams(1);
    pp.query("writeStreams",writeStreams);
    std::string streamfile = outfile+"_stream";
    pp.query("streamfile",streamfile);
    //
    int writeStreamBin(1);
    pp.query("writeStreamBin",writeStreamBin);
    std::string streamBinfile = outfile+"_streamBin";
    pp.query("streamBinfile",streamBinfile);
    
    // sanity checks
    if (inVarNames.size()!=DEF_FCOMP) 
      Abort("inVarNames wrong size");
    if (outVarNames.size()!=DEF_FCOMP) 
      Abort("outVarNames wrong size");
    // AJA - hardwired variables - inVarNames moved to StreamPC.H
    Print() << "infile = " << infile << std::endl;
    Print() << "Variables:" << std::endl;
    for (int i=0; i<inVarNames.size(); i++) 
      Print() << "  " << i << " " << inVarNames[i] << std::endl;
    Print() << "outfile = " << outfile << std::endl;
    Print() << "Variables:" << std::endl;
    for (int i=0; i<outVarNames.size(); i++) 
      Print() << "  " << i << " " << outVarNames[i] << std::endl;

    if (writeParticles)
      Print() << "Will write particles to " << partfile << std::endl;
    else
      Print() << "(Not writing particles)" << std::endl;

    if (writeStreams)
      Print() << "Will write streams to " << streamfile << std::endl;
    else
      Print() << "(Not writing streams)" << std::endl;

    if (writeStreamBin)
      Print() << "Will write binary streams to " << streamfile << std::endl;
    else
      Print() << "(Not writing binary streams)" << std::endl;

    //
    //
    //
    PlotFileData pf(infile);
    int finestLevel = pf.finestLevel();
    Vector<Geometry> geoms(finestLevel+1);
    Vector<BoxArray> grids(finestLevel+1);
    Vector<DistributionMapping> dms(finestLevel+1);
    Vector<int> ratios(finestLevel);

    Vector<int> pp_is_per(AMREX_SPACEDIM);
    pp.getarr("is_per",pp_is_per,0,AMREX_SPACEDIM);
    Array<int,AMREX_SPACEDIM> is_per = {D_DECL(pp_is_per[0], pp_is_per[1], pp_is_per[2])};

    RealBox rb(pf.probLo(),pf.probHi());

    int Nlev = finestLevel + 1;
    Vector<Vector<MultiFab>> pfdata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      geoms[lev].define(pf.probDomain(lev),rb,pf.coordSys(),is_per);
      grids[lev] = pf.boxArray(lev);
      dms[lev] = pf.DistributionMap(lev);
      if (lev < finestLevel) ratios[lev] = pf.refRatio(lev);

      pfdata[lev].resize(DEF_FCOMP);
      for (int d=0; d<DEF_FCOMP; ++d) {
        pfdata[lev][d] = pf.get(lev,inVarNames[d]);
      }
    }

    Real time=0;
    PhysBCFunctNoOp f;
    PCInterp cbi;
    BCRec bc;
    int nGrow = 3;
    pp.query("nGrow",nGrow);
    AMREX_ALWAYS_ASSERT(nGrow>=1);
    int nComp = inVarNames.size();
    Vector<MultiFab> vectorField(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      vectorField[lev].define(grids[lev],dms[lev],nComp,nGrow);
      for (int iComp=0; iComp<DEF_FCOMP; ++iComp) {
        if (lev==0) {
          FillPatchSingleLevel(vectorField[lev],time,
			       {&pfdata[lev][iComp]},{time},0,iComp,1,geoms[0],f,0);
        }
        else
        {
          FillPatchTwoLevels(vectorField[lev],time,{&pfdata[lev-1][iComp]},{time},
			     {&pfdata[lev][iComp]},{time},0,iComp,1,
                             geoms[lev-1],geoms[lev],f,0,f,0,
			     ratios[lev-1]*IntVect::Unit,&cbi,{bc},0);
        }
      }
      vectorField[lev].FillBoundary(geoms[lev].periodicity());
    }

    int Nsteps = 50;
    pp.query("Nsteps",Nsteps);
    StreamParticleContainer spc(Nsteps+1,geoms,dms,grids,ratios);

    // Get seed locations
    Print() << "Getting seed locations..." << std::endl;
    Vector<int> faceData;
    auto locs = GetSeedLocations(spc,faceData);
    int nStreamPairs = locs.size();

    // Initialise particles
    Print() << "Initiating particles..." << std::endl;
    spc.InitParticles(locs);

    // Check initialisation went ok
    spc.InspectParticles(nStreamPairs);

    // Interpolate at start
    Print() << "Interpolation at the seed points..." << std::endl;
    spc.InterpDataAtLocation(0,vectorField);

    // Check still ok
    if (!spc.OK())
      Print() << "spc not OK (before)" << std::endl;
    else 
      Print() << "spc OK (before)" << std::endl;

    // Follow streams, interpolating as we go
    Print() << "Computing streams and interpolating..." << std::endl;
    Real hRK = 0.1; pp.query("hRK",hRK);
    AMREX_ALWAYS_ASSERT(hRK>0 && hRK<=0.5);
    Real dt = hRK * geoms[finestLevel].CellSize()[0];

    for (int step=0; step<Nsteps; ++step) {

      // find next location
      spc.ComputeNextLocation(step,dt,vectorField);

      // interpolate all data
      spc.InterpDataAtLocation(step+1,vectorField);
#if 0
      // check still ok 
      if (!spc.OK())
	Print() << "spc not OK (during; step = " << step << ")" << std::endl;
#endif
    }
    
    // check in again
    spc.InspectParticles(nStreamPairs);
    if (!spc.OK())
      Print() << "spc not OK (after)" << std::endl;
    else 
      Print() << "spc OK (after)" << std::endl;
    
    //
    // Write particles
    //
    if (writeParticles) {
      spc.WritePlotFile(partfile, "particles");
    }

    //
    // Write streams
    //
    if (writeStreams) {
      Print() << "Writing streamlines in Tecplot ascii format to " << streamfile << std::endl;
      spc.WriteStreamAsTecplot(streamfile);
    }

    //
    // Write streamBin
    //
    if (writeStreamBin) {
      Print() << "Writing streamlines as binary " << streamBinfile << std::endl;
      spc.WriteStreamAsBinary(streamBinfile,faceData,nStreamPairs);
    }
    
  }
  Finalize();
  return 0;
}
