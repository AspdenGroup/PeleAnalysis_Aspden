#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "StreamPC.H"

using namespace amrex;

typedef Array<Real,AMREX_SPACEDIM> dim3;
typedef Array<int,AMREX_SPACEDIM> int3;

StreamParticleContainer::
StreamParticleContainer(int                                 a_nPtsOnStrm,
                        const Vector<Geometry>            & a_geoms,
                        const Vector<DistributionMapping> & a_dmaps,
                        const Vector<BoxArray>            & a_bas,
                        const Vector<int>                 & a_rrs)
  : ParticleContainer<0, DEF_PCOMP, 0, 0> (a_geoms, a_dmaps, a_bas, a_rrs)
{
  Nlev = a_geoms.size();
  nPtsOnStrm = a_nPtsOnStrm;
  sizeOfRealStreamData = nPtsOnStrm * DEF_PCOMP;
  for (int i=0; i<sizeOfRealStreamData; ++i) {
    AddRealComp(true);
  }
  for (int lev=0; lev<numLevels(); ++lev)
  {
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
      auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
    }
  }
}

void
StreamParticleContainer::
InitParticles (const Vector<Vector<Real>>& locs)
{
  BL_PROFILE("StreamParticleContainer::InitParticles");

  int streamLoc = 0;
  int offset = DEF_PCOMP*streamLoc;

  int lev = 0;
  int grid_id = 0;
  int tile_id = 0;
  int owner = ParticleDistributionMap(lev)[0];
  auto& particle_tile = GetParticles(0)[std::make_pair(grid_id,tile_id)];

  if (ParallelDescriptor::MyProc() == owner)
  {
    for (const auto& loc : locs)
    {
      // Keep track of pairs of lines
      Array<Long,2> ppair = {ParticleType::NextID(), ParticleType::NextID()};

      for (int i_part=0; i_part<2; i_part++)
      {
        ParticleType p;
        p.id()  = ppair[i_part];
        p.cpu() = ParallelDescriptor::MyProc();

        AMREX_D_EXPR(p.pos(0) = loc[0],
                     p.pos(1) = loc[1],
                     p.pos(2) = loc[2]);
	
        p.idata(0) = streamLoc;                  // Current position
        p.idata(1) = i_part==0 ? +1 : -1;        // Direction of integration
        p.idata(2) = ppair[ i_part==0 ? 1 : 0];  // Other line from this seed

        particle_tile.push_back(p);

        auto& soa = particle_tile.GetStructOfArrays();
        for (int i=0; i<NumRuntimeRealComps(); ++i) {
	  soa.GetRealData(i).push_back(i<AMREX_SPACEDIM ? loc[i] : 0);
	}
      }
    }
  }
  Redistribute();
}

void
StreamParticleContainer::
SetParticleLocation(int a_streamLoc, int a_nGrow)
{
  BL_PROFILE("StreamParticleContainer::SetParticleLocation");

  AMREX_ALWAYS_ASSERT(a_nGrow > 0);
  bool redist = false;
  int offset = DEF_PCOMP * a_streamLoc;
  dim3 newpos, blo, bhi;
  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx   = geom.CellSize();
    const auto& plo  = geom.ProbLo();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const auto vbx = grow(pti.validbox(),a_nGrow-1);
      const auto& vse = vbx.smallEnd();
      const auto& vbe = vbx.bigEnd();
      blo = {AMREX_D_DECL(plo[0] + vse[0]*dx[0],
                          plo[1] + vse[1]*dx[1],
                          plo[2] + vse[2]*dx[2])};

      bhi = {AMREX_D_DECL(plo[0] + (vbe[0]+1)*dx[0],
                          plo[1] + (vbe[1]+1)*dx[1],
                          plo[2] + (vbe[2]+1)*dx[2])};

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
        if (p.id() > 0)
        {
	  int pId = p.id();

          newpos = {AMREX_D_DECL(soa.GetRealData(offset + RealData::xloc)[pindex],
                                 soa.GetRealData(offset + RealData::yloc)[pindex],
                                 soa.GetRealData(offset + RealData::zloc)[pindex])};

          AMREX_D_EXPR(p.pos(0) = newpos[0], p.pos(1) = newpos[1], p.pos(2) = newpos[2]);

          for (int d=0; d<AMREX_SPACEDIM; ++d)
          {
            redist |= (newpos[d]<blo[d] || newpos[d]>bhi[d]);
          }
        }
      }
    }
  }

  ParallelDescriptor::ReduceBoolOr(redist);
  if (redist) {
    Redistribute();
  }
}

static void vnrml(Vector<Real>& vec, int dir)
{
  static Real eps = 1.e-12;
  Real sum = AMREX_D_TERM(  vec[0] * vec[0],
                          + vec[1] * vec[1],
                          + vec[2] * vec[2]);
  Vector<Real> u = vec;
  if (sum > eps) {
    sum = 1. / std::sqrt(sum);
    for (int i=0; i<AMREX_SPACEDIM; ++i) vec[i] *= dir * sum;
  }
  else {
    vec = {AMREX_D_DECL(0, 0, 0)};
  }
}

static bool ntrpv(const dim3& x,const FArrayBox& gfab,
                  const Real* dx,const Real* plo,const Real* phi,Vector<Real>& u,int nComp)
{
  int3 b; // node based index of lower-left corner of cell
  dim3 n; // interpolation factor in [0,1]^3

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    b[d] = std::floor( (x[d] - plo[d]) / dx[d] - 0.5 );
    n[d] = ( x[d] - ( (b[d] + 0.5 ) * dx[d] + plo[d] ) )/dx[d];
    n[d] = std::max(0., std::min(1.,n[d]));
  }

  const auto& gbx = gfab.box();
  const auto& glo = gbx.smallEnd();
  const auto& ghi = gbx.bigEnd();
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (b[d] < glo[d] ||  b[d] > ghi[d]-1) {
      Print() << "dir: " << d << std::endl;
      Print() << "d,b,glo,ghi: " << d << " " << b[d] << " " << glo[d] << " " << ghi[d] << std::endl;
      Print() << "x,plo,phi " << x[d] << " " << plo[d] << " " << phi[d] << std::endl;
      Print() << "boxlo,boxhi " << plo[d]+glo[d]*dx[d] << " " << plo[d]+(ghi[d]+1)*dx[d] << std::endl;
      return false;
    }
  }

  const auto& g = gfab.array();
  for (int i=0; i<nComp; ++i) {
#if AMREX_SPACEDIM == 2
    u[i] =
      +   n[0]   *   n[1]   * g(b[0]+1,b[1]+1,0,i)
      +   n[0]   * (1-n[1]) * g(b[0]+1,b[1]  ,0,i)
      + (1-n[0]) *   n[1]   * g(b[0]  ,b[1]+1,0,i)
      + (1-n[0]) * (1-n[1]) * g(b[0]  ,b[1]  ,0,i);
#else
    u[i] =
      +    n[0]   *    n[1]  *    n[2]  * g(b[0]+1,b[1]+1,b[2]+1,i)
      +    n[0]   * (1-n[1]) *    n[2]  * g(b[0]+1,b[1]  ,b[2]+1,i)
      +    n[0]   *    n[1]  * (1-n[2]) * g(b[0]+1,b[1]+1,b[2]  ,i)
      +    n[0]   * (1-n[1]) * (1-n[2]) * g(b[0]+1,b[1]  ,b[2]  ,i)
      +  (1-n[0]) *    n[1]  *    n[2]  * g(b[0]  ,b[1]+1,b[2]+1,i)
      +  (1-n[0]) * (1-n[1]) *    n[2]  * g(b[0]  ,b[1]  ,b[2]+1,i)
      +  (1-n[0]) *    n[1]  * (1-n[2]) * g(b[0]  ,b[1]+1,b[2]  ,i)
      +  (1-n[0]) * (1-n[1]) * (1-n[2]) * g(b[0]  ,b[1]  ,b[2]  ,i);
#endif
  }
  return true;
}

static bool
RK4(dim3 & x,Real dt,const FArrayBox& v,const Real* dx,const Real* plo,const Real* phi,int dir)
{
  Vector<Real> vec(AMREX_SPACEDIM);
  dim3 k1, k2, k3, k4;
  dim3 xx = x;
  
  if (!ntrpv(xx,v,dx,plo,phi,vec,AMREX_SPACEDIM)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k1[d] = vec[d] * dt;
    xx[d] = x[d] + k1[d] * 0.5;
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec,AMREX_SPACEDIM)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k2[d] = vec[d] * dt;
    xx[d] = x[d] + k2[d] * 0.5;
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec,AMREX_SPACEDIM)) return false;
  vnrml(vec,dir);

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k3[d] = vec[d] * dt;
    xx[d] = x[d] + k3[d];
  }
  if (!ntrpv(xx,v,dx,plo,phi,vec,AMREX_SPACEDIM)) return false;
  vnrml(vec,dir);

  const Real third = 1./3.;
  const Real sixth = 1./6.;
  dim3 delta;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    k4[d] = vec[d] * dt;
    delta[d] = (k1[d] + k4[d])*sixth + (k2[d] + k3[d])*third;
  }

  // just do basic update without dealing with domain...
  for (int d=0; d<AMREX_SPACEDIM; ++d)
    x[d] += delta[d];
  
  /* Old way to control domain; hacked out by aja

  // cut step length to keep in domain (FIXME: Deal with periodic)
  Real scale = 1;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (x[d]+delta[d] < plo[d]) {
      scale = std::min(scale, std::abs((x[d] - plo[d])/delta[d]));
    }
    if (x[d]+delta[d] > plo[d]) {
      scale = std::min(scale, std::abs((phi[d] - x[d])/delta[d]));
    }
  }
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    x[d] += scale * delta[d];
    x[d] = std::min(phi[d]-1.e-10, std::max(plo[d]+1.e-10, x[d]) ); // Deal with precision issues
  }

  */
  
  return true;
}

void
StreamParticleContainer::
ComputeNextLocation(int                      a_fromLoc,
                    Real                     a_delta_t,
                    const Vector<MultiFab> & a_vectorField)
{
  BL_PROFILE("StreamParticleContainer::ComputeNextLocation");

  const int nGrow = a_vectorField[0].nGrow();
  SetParticleLocation(a_fromLoc,nGrow);

  const int new_loc_id = a_fromLoc + 1;
  int offset = DEF_PCOMP * new_loc_id;

  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& v = a_vectorField[lev][pti];

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
	
	int pId = p.id();
	
	const int dir = p.idata(1);
        dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
        if (p.id() > 0)
        {
          if (!RK4(x,a_delta_t,v,dx,plo,phi,dir))
          {
            Abort("bad RK");
          }
        }
	// put particle position into stream
        AMREX_D_EXPR(soa.GetRealData(offset + RealData::xloc)[pindex] = x[0],
                     soa.GetRealData(offset + RealData::yloc)[pindex] = x[1],
                     soa.GetRealData(offset + RealData::zloc)[pindex] = x[2]);
      }
    }
  }
}

void
StreamParticleContainer::
InterpDataAtLocation(int                      a_fromLoc,
		     const Vector<MultiFab> & a_vectorField)
{
  BL_PROFILE("StreamParticleContainer::InterpDataAtLocation");

  const int nGrow = a_vectorField[0].nGrow();
  SetParticleLocation(a_fromLoc,nGrow);

  int offset = DEF_PCOMP * a_fromLoc; // components on particle

  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    dim3 Lx;
    for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++)
      Lx[iComp] = phi[iComp]-plo[iComp];

    //Print() << "Lx = " << Lx << std::endl;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& v = a_vectorField[lev][pti];

      for (int pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
	if (p.id()>0) {
	  int pId = p.id();

	  // where's the particle?
	  dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
	  Vector<Real> ntrpvOut(DEF_FCOMP); // components in infile

	  // interpolate all data to particle location
	  ntrpv(x,v,dx,plo,phi,ntrpvOut,DEF_FCOMP); // components in infile

	  // copy the interpolated data to the particle
	  // first DIM components are particle location
	  // next DIM components are stream location w/o adjusting for periodicity
	  // then we have the interpolated data we want
	  for (int iComp=AMREX_SPACEDIM; iComp<DEF_FCOMP; ++iComp) {
	    int idxOnPart = offset + iComp + AMREX_SPACEDIM;
	    soa.GetRealData(idxOnPart)[pindex] = ntrpvOut[iComp];
	  }
	  // let's figure out the location on the stream w/o adjusting for periodicity
	  if (a_fromLoc==0) { // nothing to do on the surface; just take a copy of the location
	    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	      int idx = offset + d;
	      soa.GetRealData(idx+AMREX_SPACEDIM)[pindex] = soa.GetRealData(idx)[pindex];
	    }
	  } else { // set new location adjusting for periodicity
	    // calculate the change in position delta
	    // add to the old position
	    // adjust delta if it's affected by periodicity (i.e. delta too big)
	    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	      int  idx    = offset + d;
	      int  idxOld = idx - DEF_PCOMP;
	      Real xnew   = soa.GetRealData(idx   )[pindex];
	      Real xold   = soa.GetRealData(idxOld)[pindex];
	      Real delta  = xnew-xold;
	      if (fabs(delta)>dx[d]) { // has been adjusted for periodicity
		//printf("%i %i %e %e %e",pindex,d,xnew,xold,delta);
		if (delta<0.) delta+=Lx[d];
		else          delta-=Lx[d];
		//printf(" --> %e\n",delta);
	      }
	      // store periodicity-adjusted copy of location at idx+SPACEDIM
	      Real sold   = soa.GetRealData(idxOld+AMREX_SPACEDIM)[pindex];
	      Real snew   = sold+delta;
	      soa.GetRealData(idx+AMREX_SPACEDIM)[pindex] = snew;
	    }
	  }
	  //hack last component to AMR level for debugging
	  //soa.GetRealData(offset + DEF_FCOMP-1)[pindex] = (Real)lev;
	}
      }
    }
  }
}

void
StreamParticleContainer::
WriteStreamAsTecplot(const std::string& outFile)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0,1);

  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (int pindex=0; pindex<aos.size() && !will_write; ++pindex)
      {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }
    }
  }

  if (will_write)
  {
    std::string fileName = outFile + "/str_";
    fileName = Concatenate(fileName,myProc) + ".dat";
    std::ofstream ofs(fileName.c_str());

    ofs << "VARIABLES = ";
    for (int iComp=0; iComp<DEF_FCOMP; ++iComp)
      ofs << outVarNames[iComp] << " ";
    ofs << '\n';
    
    int minId=100000000;
    int maxId=-minId;
    for (int lev = 0; lev < Nlev; ++lev)
    {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();

        for (int pindex=0; pindex<aos.size(); ++pindex)
        {
	  ParticleType& p = aos[pindex];
	  if (p.id()>0) {
	    int pId = p.id();
	    minId=min(minId,pId);
	    maxId=max(maxId,pId);

	    ofs << "ZONE I=1 J=" << nPtsOnStrm << " K=1 FORMAT=POINT\n";
	    for (int j=0; j<nPtsOnStrm; ++j)
	      {
		// by including the spacedim offset, we use locations w/o periodicity adjustments
		int offset = j*DEF_PCOMP + AMREX_SPACEDIM;
		for (int iComp=0; iComp<DEF_FCOMP; ++iComp) {
		  ofs << soa.GetRealData(offset + iComp)[pindex] << " ";
		}
		ofs << '\n';
	      }
	  }
        }
      }
    }
    ofs.close();
    ParallelDescriptor::ReduceIntMin(minId);
    ParallelDescriptor::ReduceIntMax(maxId);
    //if (ParallelDescriptor::IOProcessor())
    //printf("writeTec: minId / maxId = %i / %i\n",minId,maxId);
  }
}

void
StreamParticleContainer::
WriteStreamAsBinary(const std::string& outFile, Vector<int>& faceData)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0,1);

  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (int pindex=0; pindex<aos.size() && !will_write; ++pindex)
      {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }

    }
  }
  
  // Need to count the total number of streams to be written
  // by all ptiters on all levels on this processor
  int nStreams = 0;
  int nStreamsCheck = 0;
  
  if (will_write)
  {
    for (int lev = 0; lev < Nlev; ++lev) {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
	auto& aos = pti.GetArrayOfStructs();
	for (int pindex=0; pindex<aos.size(); ++pindex) {
	  ParticleType& p = aos[pindex];
	  if (p.id() > 0) {
	    nStreams++;
	  }
	}
      }
    }
    
    // write to a binary file
    std::string rootName = outFile + "/str_";
    std::string fileName = Concatenate(rootName,myProc) + ".bin";
    std::string headName = Concatenate(rootName,myProc) + ".head";
    FILE *file=fopen(fileName.c_str(),"w");
    FILE *head=fopen(headName.c_str(),"w");
    // total number of streams in file
    fwrite(&(nStreams),sizeof(int),1,file);

    int minId=100000000;
    int maxId=-minId;
    for (int lev = 0; lev < Nlev; ++lev)
    {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();
	
        for (int pindex=0; pindex<aos.size(); ++pindex)
        {
	  ParticleType& p = aos[pindex];
	  if (p.id()>0) {
	    // write info about this stream
	    int pId = p.id();
	    fwrite(&(pId),sizeof(int),1,file);    // id
	    int dir = p.idata(1);
	    fwrite(&(dir),sizeof(int),1,file);    // dir
	    int pairId = p.idata(2);
	    fwrite(&(pairId),sizeof(int),1,file); // pair id
	    
	    fprintf(head,"%i %i %i\n",pId,dir,pairId);
	    
	    minId=min(minId,pId);
	    maxId=max(maxId,pId);
	    
	    for (int j=0; j<nPtsOnStrm; ++j) {
	      // by including the spacedim offset, we use locations w/o periodicity adjustments
	      int offset = j*DEF_FCOMP + AMREX_SPACEDIM;
	      fwrite(&(soa.GetRealData(offset)[0]),sizeof(Real),DEF_FCOMP,file);
	    }
	    
	    nStreamsCheck++; // sanity check to make sure we wrote number of streams anticipated
	  }
        }
      }
    }
    ParallelDescriptor::ReduceIntMin(minId);
    ParallelDescriptor::ReduceIntMax(maxId);
    //if (ParallelDescriptor::IOProcessor())
    //printf("writeBin: minId / maxId = %i / %i\n",minId,maxId);
    
    //std::cout << nStreams << " ?= "  << nStreamsCheck << std::endl;
    
    if (nStreams!=nStreamsCheck)
      std::cout << "(nStreams!=nStreamsCheck) : "
		<< nStreams << " != "  << nStreamsCheck << std::endl;
    
    fclose(file);
    fclose(head);
  }

  // Write header file with everything consistent across all processors
  ParallelDescriptor::ReduceIntSum(nStreams);
  if (ParallelDescriptor::IOProcessor()) {
    std::string fileName = outFile + "/Header";
    std::ofstream ofs(fileName.c_str());
    ofs << "Even odder-ball replacement for sampled streams" << std::endl;
    ofs << nProcs << std::endl;     // translates to number of files to read
    ofs << nStreams << std::endl;   // total number of streams
    ofs << nPtsOnStrm << std::endl; // number of points
    ofs << DEF_FCOMP << std::endl;  // number of variables
    for (int iComp=0; iComp<DEF_FCOMP; ++iComp)
      ofs << outVarNames[iComp] << " ";
    ofs << std::endl;
    ofs << faceData.size() << std::endl;
    ofs.write((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ofs << '\n';
    ofs.close();
  }

}


void
StreamParticleContainer::
InspectParticles (const int nStreamPairs)
{
  BL_PROFILE("StreamParticleContainer::InspectParticles");

  SetParticleLocation(0,1);
  
  int IOProc = ParallelDescriptor::IOProcessor();
  int myProc = ParallelDescriptor::MyProc();
  int nProcs = ParallelDescriptor::NProcs();

  int partIndexing[nProcs][2*nStreamPairs+1];
  int ptiIndexing[nProcs][2*nStreamPairs+1];

  for (int iProc = 0; iProc<nProcs; iProc++) {
    for (int iStream=0; iStream<=2*nStreamPairs; iStream++) {
      partIndexing[iProc][iStream] = -1;
      ptiIndexing[iProc][iStream]  = -1;
    }
  }

  int ptiCounter;
  
  int minId=100000000;
  int maxId=-minId;

  for (int lev = 0; lev < Nlev; ++lev) {
    ptiCounter=0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti, ++ptiCounter) {

      auto& aos = pti.GetArrayOfStructs();

      int NP = aos.size();

      for (int pindex=0; pindex<aos.size(); ++pindex) {
        ParticleType& p = aos[pindex];

        if (p.id() > 0) {
	  int pId = p.id();
	  minId=min(minId,pId);
	  maxId=max(maxId,pId);
	  if (pId>2*nStreamPairs) {
	    std::cout << "pId > 2*nStreamPairs " << pId << " > "
		      << (2*nStreamPairs) << std::endl;
	  }
	  partIndexing[myProc][pId] = pindex;
	  ptiIndexing[myProc][pId] = ptiCounter;
	} else {
	  Print() << "p.id() = " << p.id() << " !> 0" <<std::endl;
	} // p.id

      } // pindex
    } // pti
    
    ParallelDescriptor::Barrier();
    // Now let's see if the pair is on the same processor and pti
    for (int iProc=0; iProc < nProcs; iProc++) {
      if (iProc==myProc) {
	ptiCounter=0; // reset
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti, ++ptiCounter) {
	  auto& aos = pti.GetArrayOfStructs();
	  int goodCount=0;
	  for (int iStream=1; iStream<=2*nStreamPairs; iStream++) {
	    int pindex1 = partIndexing[myProc][iStream];
	    int ptiindex1 = ptiIndexing[myProc][iStream];
	    //if ( (pindex1>-1) && (pindex1<aos.size()) ) {
	    if ( (pindex1>-1) && (ptiindex1==ptiCounter) ) {
	      ParticleType& p1 = aos[pindex1];
	      int myId1   = p1.id();
	      int pairId1 = p1.idata(2);
	      int pindex2 = partIndexing[myProc][pairId1];
	      ParticleType& p2 = aos[pindex2];
	      int myId2   = p2.id();
	      int pairId2 = p2.idata(2);
	      if ( (myId1!=pairId2) || (pairId1!=myId2) ) {
		std::cout << "myProc / pindex1 / myId1 / pairId1 = " 
			  << myProc << " / "
			  << pindex1 << " / "
			  << myId1 << " / "
			  << pairId1 << " / " << std::endl;
		std::cout << "myProc / pindex2 / myId2 / pairId2 = "
			  << myProc << " / "
			  << pindex2 << " / "
			  << myId2 << " / "
			  << pairId2 << " / " << std::endl;
	      } else {
		// looks good
		goodCount++;
	      }
	    }
	  }
	  if (goodCount!=aos.size()) {
	    std::cout << "aos.size() != goodCount : "
		      << aos.size() << " != " 
		      << goodCount << std::endl;
	  }
	}
      }
      ParallelDescriptor::Barrier();
    }
    
  } // lev

  ParallelDescriptor::ReduceIntMin(minId);
  ParallelDescriptor::ReduceIntMax(maxId);
  if (ParallelDescriptor::IOProcessor())
    printf("InspectParticles: minId / maxId = %i / %i\n",minId,maxId);
  
  int duplicate(0);
  int missing(0);
  for (int iStream=1; iStream<=2*nStreamPairs; iStream++) {
    int pIdx   = partIndexing[myProc][iStream];
    int ptiIdx = ptiIndexing[myProc][iStream];
    int gotcha=0;
    if ( (pIdx!=-1) && (ptiIdx!=-1) ) {
      gotcha++;
    }
    ParallelDescriptor::ReduceIntSum(gotcha);
    if (gotcha==0) {
      Print() << "iStream(" << iStream << ") == 0" << std::endl;
      missing++;
    }
    if (gotcha>1) {
      Print() << "iStream(" << iStream << ") > 1" << std::endl;
      duplicate++;
    }
  }
  if (missing>0)
    Print() << "Missing " << missing << " particles !!!" << std::endl;
  else {
    Print() << "No particles missing :-)" << std::endl;
  }    
  if (duplicate>0)
    Print() << "There are " << duplicate << " duplicates !!!" << std::endl;
  else {
    Print() << "No duplicate particles :-)" << std::endl;
  }
  
}
