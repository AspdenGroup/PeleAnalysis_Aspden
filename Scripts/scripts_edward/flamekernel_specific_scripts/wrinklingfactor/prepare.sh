#!/bin/bash

export dim=2
export loc=~/D1_SSD/AMReX/PeleAnalysis_Aspden;

ln -s ${loc}/Src/plotProg${dim}d.gnu.MPI.ex .;
ln -s ${loc}/Src/isosurface${dim}d.gnu.MPI.ex .;
ln -s ${loc}/Src/areaMEF${dim}d.gnu.ex .;
ln -s ${loc}/AmrDerive/AmrDeriveProgIntegration${dim}d.Linux.g++.gfortran.ex .;
