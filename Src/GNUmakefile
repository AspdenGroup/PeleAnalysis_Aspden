AMREX_HOME = ../amrex
DEBUG	      = FALSE
DIM	      = 3
COMP          = gnu
PRECISION     = DOUBLE
USE_MPI       = FALSE
USE_OMP       = FALSE

EBASE         = combinePlts
EBASE         = isosurface
EBASE         = partStream
#EBASE         = stream2plt
#EBASE         = streamScatter
#EBASE         = streamSub
#EBASE         = streamBinTubeStats
#EBASE         = streamBinTubeStats_TLH
#EBASE         = curvature
#EBASE         = surfMEFtoDATnorm
#EBASE         = surfMEFtoDAT
#EBASE         = surfMEFtoDATbasic
#EBASE         = surfDATtoMEF
#EBASE         = buildPMF
#EBASE         = subPlt
#EBASE         = conditionalMean
#EBASE         = jpdf
#EBASE         = plotProg
#EBASE         = plotPhi
#EBASE         = decimate
#EBASE         = areaMEF   
#EBASE		= deltaCheck
#EBASE		= avgPlots
#EBASE		= surfStats
#EBASE		= makePlotfile
#EBASE		= addDummies
#EBASE		= diffSlice
#EBASE		= centerline
#EBASE		= plotZC_noTerms
#EBASE		= flameindex
#EBASE		= trimPlotfile
#EBASE		= jpdf
#EBASE		= machnumber
#EBASE		= dissipationrateZ
#EBASE		= slicescatter
#EBASE		= avgFlucs
#EBASE		= avgPlots
#EBASE		= plotZC_noTerms
#EBASE		= boxTagging
#EBASE		= straintensor
#EBASE		= microscales
#EBASE		= diffusivefluxes
EBASE		= grad
NEEDS_f90_SRC = FALSE

#EBASE         = sampleStreamlines
#EBASE         = stream
#EBASE		= plotZC
#EBASE		= avgTheta
#NEEDS_f90_SRC = TRUE


include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif

CEXE_headers += StreamPC.H
ifeq ($(EBASE),partStream)
  CEXE_sources += StreamPC.cpp
endif

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

Pdirs   := Base Boundary AmrCore Extern/amrdata LinearSolvers/MLMG Particle
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

#CEXE_sources += AppendToPlotFile.cpp WritePlotFile.cpp
#CEXE_headers += AppendToPlotFile.H   WritePlotFile.H
INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
VPATH_LOCATIONS += $(AMREX_HOME)/Tools/C_util

CEXE_sources += AMReX_Extrapolater.cpp
CEXE_headers += AMReX_Extrapolater.H
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
VPATH_LOCATIONS += $(AMREX_HOME)/Src/Amr

SDF_LOC = ../Tools/SDFGen
include $(SDF_LOC)/Make.package
INCLUDE_LOCATIONS += $(SDF_LOC)
VPATH_LOCATIONS += $(SDF_LOC)

vpath %.c   : $(VPATH_LOCATIONS)
vpath %.h   : $(VPATH_LOCATIONS)
vpath %.cpp : $(VPATH_LOCATIONS)
vpath %.H   : $(VPATH_LOCATIONS)
vpath %.F   : $(VPATH_LOCATIONS)
vpath %.f   : $(VPATH_LOCATIONS)
vpath %.f90 : $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
