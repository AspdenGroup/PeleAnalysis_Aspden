#!/bin/bash

export isoval=0.9
export fuel=H2
export dim=2
export ntasks=2

rm -r *plt?????_prog* area.dat *~ Backtrace* int_prog.dat wrinkling_factor.csv;

for i in plt?????;
do
    ./plotProg${dim}d.gnu.MPI.ex -n ${ntasks} infile= $i;
    ./isosurface${dim}d.gnu.MPI.ex -n ${ntasks} inputs \
		comps=1 \
		isoCompName=prog_${fuel} \
		isoVal=${isoval} \
		infile=${i}_prog
done
./areaMEF${dim}d.gnu.ex inputs infile= *.mef;
./AmrDeriveProgIntegration${dim}d.Linux.g++.gfortran.ex infile= plt?????_prog vars=prog_${fuel};
python wf.py;

rm -r plt?????_prog* area.dat int_prog.dat;
