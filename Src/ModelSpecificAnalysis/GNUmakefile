TOP = .
PELE_ANALYSIS_HOME = ../..
PELE_PHYSICS_HOME = PelePhysics
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
EBASE = plotTransportCoeff

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
#Eos_dir     := Soave-Redlich-Kwong
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego
#Reactions_dir := Null

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
#Chemistry_Model := Davis
Chemistry_Model := LiDryer
#Chemistry_Model := grimech12
#Chemistry_Model := drm19

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
#Transport_dir := EGLib
Transport_dir := Simple

include Make.ModelSpecificAnalysis
