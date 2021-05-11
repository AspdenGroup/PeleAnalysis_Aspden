AMREX_HOME = ../../../amrex

DEBUG	      = FALSE
DIM	      = 3
COMP          = gcc
PRECISION     = DOUBLE
USE_MPI       = FALSE
USE_OMP       = FALSE
EBASE         = qslim

MIX_ROOT       = ../qslim-2.1/mixkit/src
GFX_ROOT       = ../qslim-2.1/libgfx
COMP           = gcc

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DAJASUPPRESS

CEXE_sources += cmdline.cpp main.cpp output.cpp qslim.cpp

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

Pdirs   := Base Boundary AmrCore Extern/amrdata
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

CEXE_sources += cmdline.cpp main.cpp output.cpp qslim.cpp
INCLUDE_LOCATIONS += ${MIX_ROOT}
INCLUDE_LOCATIONS += ${GFX_ROOT}/include
LIBRARY_LOCATIONS += ${MIX_ROOT}
LIBRARY_LOCATIONS += ${GFX_ROOT}/src
LIBRARIES += -lmix -lgfx 
DEFINES += -DHAVE_CONFIG_H -DHAVE_BOOL

vpath %.c   : $(VPATH_LOCATIONS)
vpath %.h   : $(VPATH_LOCATIONS)
vpath %.cpp : $(VPATH_LOCATIONS)
vpath %.H   : $(VPATH_LOCATIONS)
vpath %.F   : $(VPATH_LOCATIONS)
vpath %.f   : $(VPATH_LOCATIONS)
vpath %.f90 : $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

