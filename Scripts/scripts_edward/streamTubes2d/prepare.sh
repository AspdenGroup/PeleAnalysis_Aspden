#!/bin/bash

export dim=2;
export loc=~/src/PeleAnalysis_Aspden/Src;

ln -s ${loc}/plotProg${dim}d.gnu.MPI.ex .;
ln -s ${loc}/grad${dim}d.gnu.MPI.ex .;
ln -s ${loc}/combinePlts${dim}d.gnu.MPI.ex .;
ln -s ${loc}/curvature${dim}d.gnu.MPI.ex .;
ln -s ${loc}/isosurface${dim}d.gnu.MPI.ex .;
ln -s ${loc}/stream${dim}d.gnu.MPI.ex .;
ln -s ${loc}/sampleStreamlines${dim}d.gnu.MPI.ex .;
ln -s ${loc}/streamTubeStats${dim}d.gnu.ex .;
ln -s ${loc}/surfMEFtoDAT${dim}d.gnu.MPI.ex .;
ln -s ${loc}/surfMEFtoDATbasic${dim}d.gnu.MPI.ex .;
