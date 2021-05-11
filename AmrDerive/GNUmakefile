BOXLIB_HOME = ../BoxLib

BL_NOFAST=TRUE

#
# Variables for the user to set ...
#
########## define EBASE and BL_HAS_FORT here

#
# These don't need fortran (set BL_HAS_FORT=FALSE)
#
#EBASE = AmrDeriveBinaryDump
#EBASE = AmrDeriveRmsVel
#EBASE = AmrDeriveMakePPM

#BL_HAS_FORT = FALSE

#
# These need fortran (set BL_HAS_FORT=TRUE)
#
#EBASE = AmrDeriveFilter
#EBASE = AmrDeriveFilter2D
#EBASE = AmrDeriveGenericIntegration
EBASE = AmrDeriveProgIntegration
#EBASE = AmrDerive_udot_grad2u
#EBASE = AmrDeriveCoarsen
#EBASE = AmrDerive_udot_F
#EBASE = AmrDeriveHackVelocity
#EBASE = AmrDeriveBinaryDump
#EBASE = AmrDerive_FlameDensity
#EBASE = AmrDeriveH2OIntegration

BL_HAS_FORT = TRUE

#EBASE = AmrDeriveSpectrum
#EBASE  = AmrDeriveMakePPM
#BL_HAS_FORT = FALSE
#
#
#

DEBUG	  = FALSE
PROFILE   = FALSE
DIM       = 2
USE_MPI   = FALSE
COMP      = g++
FCOMP     = gfortran
BUILDCHEM = FALSE

USE_THREADS   = FALSE
USE_ARRAYVIEW = FALSE
AIXLARGEMEM   = FALSE

HERE = .

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

DEFINES           += -DBL_NOLINEVALUES

CEXE_sources += $(EBASE).cpp
ifeq ($(EBASE), AmrDeriveQuench)
  CHEMDIR = $(BOXLIB_HOME)/HeatTransfer/ChemKinDriver
endif

include ./Make.package

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib

include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
INCLUDES += $(BOXLIB_HOME)/Src/C_BoundaryLib

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_AMRLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata
PATHDIRS += $(BOXLIB_HOME)/Src/C_BoundaryLib

ifeq (${BL_HAS_FORT}, TRUE)
FEXE_sources += ${EBASE}_F.F
endif

ifeq ($(EBASE), AmrDeriveFilter)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F FILCC_2D.F
endif

ifeq ($(EBASE), AmrDeriveSpectrum)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 LIBRARIES += -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
endif

ifeq ($(EBASE), AmrDeriveGenericIntegration)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F
 ifeq ($(DIM), 2)
  FEXE_sources += FILCC_2D.F
 endif
 ifeq ($(DIM), 3)
  FEXE_sources += FILCC_3D.F
 endif
endif

ifeq ($(EBASE), AmrDeriveHackVelocity)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F
 ifeq ($(DIM), 2)
  FEXE_sources += FILCC_2D.F
 endif
 ifeq ($(DIM), 3)
  FEXE_sources += FILCC_3D.F
 endif
endif



ifeq ($(EBASE), AmrDerive_FlameDensity)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F
 ifeq ($(DIM),2)
  FEXE_sources += FILCC_2D.F
 endif
 ifeq ($(DIM),3)
  FEXE_sources += FILCC_3D.F
 endif
endif





ifeq ($(EBASE), AmrDeriveHit)
ifeq ($(USE_MPI), TRUE)
 LIBRARIES += -lfftw_mpi -lfftw
 DEFINES += -DDO_SPECTRA_PARALLEL
else
 LIBRARIES += -lfftw3
 DEFINES += -DDO_SPECTRA_SERIAL
endif
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F FILCC_3D.F
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

$(objEXETempDir)/AmrDeriveIso_F.o: AmrDeriveIso_F.F
	@if [ ! -d $(objEXETempDir) ]; then mkdir -p $(objEXETempDir); fi
	@if [ ! -d $(f77TempDir) ]; then mkdir -p $(f77TempDir); fi
	$(FORT_CPP) -DBL_LANG_FORT $(CPPFLAGS) $< | $(FORTPREP) > $(f77TempDir)/AmrDeriveIso_F.f
	$(FC) -O0 $(fincludes) -c $(FEXEFLAGS) $(f77TempDir)/AmrDeriveIso_F.f  -o $(objEXETempDir)/AmrDeriveIso_F.o

$(objEXETempDir)/AmrDeriveRmsVel.o: AmrDeriveRmsVel.cpp
	@echo Compiling $*.cpp ...
	@if [ ! -d $(objEXETempDir) ]; then mkdir -p $(objEXETempDir); fi
	$(SILENT) $(CCACHE) $(CXX) $(CXXFLAGS) -O0 $(CPPFLAGS) -c $(CXXEXEFLAGS) AmrDeriveRmsVel.cpp -o o/3d.Linux.g++.gfortran.EXE/AmrDeriveRmsVel.o

$(objEXETempDir)/AmrDeriveMakePPM.o: AmrDeriveMakePPM.cpp
	@echo Compiling $*.cpp ...
	@if [ ! -d $(objEXETempDir) ]; then mkdir -p $(objEXETempDir); fi
	$(SILENT) $(CCACHE) $(CXX) $(CXXFLAGS) -O0 $(CPPFLAGS) -c $(CXXEXEFLAGS) AmrDeriveMakePPM.cpp -o $(objEXETempDir)/AmrDeriveMakePPM.o

ifeq ($(EBASE), AmrDerive_udot_F)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
 include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
 FEXE_sources += xtra_F.F FILCC_3D.F
endif
