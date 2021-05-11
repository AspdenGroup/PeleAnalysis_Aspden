TOP = .
PELE_ANALYSIS_HOME = ../..
PELE_PHYSICS_HOME = ../../../PeleProduction_ISC2022/Submodules/PelePhysics
PELELM_HOME = ../../../PeleProduction_ISC2022/Submodules/PeleLM
AMREX_HOME = ../../amrex
PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

#DIM        = 1
#DIM        = 2
DIM        = 3

COMP	   = gcc
USE_MPI    = TRUE
USE_OMP    = FALSE

EBASE = plotYtoX
EBASE = plotXtoY
EBASE = testQPDtools
EBASE = plotQPD
EBASE = plotProg

USE_SUNDIALS_PP = FALSE

UTIL_LOC = ../../Tools/Util
include $(UTIL_LOC)/Make.package
include $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Make.package
# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
#Eos_dir     := Soave-Redlich-Kwong
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego
#Reactions_dir := Null

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
#Chemistry_Model := Davis
Chemistry_Model := BurkeDryer_mod_noArHe
#Chemistry_Model := grimech12
#Chemistry_Model := drm19

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Simple

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
