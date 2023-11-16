#ifndef STREAM_PARTICLE_CONTAINER_H_
#define STREAM_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

// file variables in/out
#define DEF_FCOMP (AMREX_SPACEDIM+7)

// variables to carry on the stream
// keep an extra copy of location to deal with periodic bcs
#define DEF_PCOMP (AMREX_SPACEDIM+DEF_FCOMP)

const amrex::Vector<std::string> inVarNames = {
  //  D_DECL("x_velocity", "y_velocity", "z_velocity"),
  //D_DECL("FlameNormalX_prog_H2", "FlameNormalY_prog_H2", "FlameNormalZ_prog_H2"),
  /*"Y(H2)",
  "temp",
  "H2_ConsumptionRate",
  "HeatRelease",
  "StrainRate_prog_H2",
  "||gradtemp||",
  "MeanCurvature_prog_H2",
  "GaussianCurvature_prog_H2",
  "phi"*/
  D_DECL("FlameNormalX_progress_variable", "FlameNormalY_progress_variable", "FlameNormalZ_progress_variable"),  
  "Y(H2)",
  "temp",
  "density",
  "I_R(H2)",
  "QoverCP",
  "phi",
  "zone"
}; 

const amrex::Vector<std::string> outVarNames = {
  D_DECL("X", "Y", "Z"),
  "H2_MassFraction",
  "Temp",
  "Density",
  "H2_ProductionRate",
  "QoverCP",
  "phi",
  "zone"
}; 

struct RealData
{
  enum {
    AMREX_D_DECL(xloc = 0,
                 yloc,
                 zloc)
  };
};

class StreamParticleContainer
  : public amrex::ParticleContainer<0, DEF_PCOMP, 0, 0>
{
public:

  using MyParIter = amrex::ParIter<0, DEF_PCOMP, 0, 0>;

  StreamParticleContainer (int                                               a_nPtsOnStrm,
                           const amrex::Vector<amrex::Geometry>            & a_geoms,
                           const amrex::Vector<amrex::DistributionMapping> & a_dmaps,
                           const amrex::Vector<amrex::BoxArray>            & a_bas,
                           const amrex::Vector<int>                        & a_rrs);

  void InitParticles(const amrex::Vector<amrex::Vector<amrex::Real>>& a_locs);

  void SetParticleLocation(int a_streamPos, int a_nGrow);

  void ComputeNextLocation(int                                    a_fromLoc,
                           amrex::Real                            a_delta_t,
                           const amrex::Vector<amrex::MultiFab> & a_vectorField);

  void InterpDataAtLocation(int                                    a_fromLoc,
                            const amrex::Vector<amrex::MultiFab> & a_vectorField);

  void WriteStreamAsTecplot(const std::string& outfile);

  void WriteStreamAsBinary(const std::string& outfile,
			   amrex::Vector<int>& faceData,
			   const int nStreamPairs);

  void InspectParticles(const int nStreamPairs);

protected:
  int Nlev;
  int nPtsOnStrm;
  int sizeOfRealStreamData;

private:

};

#endif
