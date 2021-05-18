#!/bin/bash

export isoval=0.9
export fuel=H2
export ntasks=2

rm -r int_prog.dat radius.csv Backtrace.* *~ *_prog;
for i in plt*;
do
    ./plotProg2d.gnu.MPI.ex -n ${ntasks} infile=${i} fuelName=${fuel};
done
./AmrDeriveProgIntegration2d.Linux.g++.gfortran.ex infile= plt?????_prog vars=prog_${fuel};
python r.py;
rm -r *_prog int_prog.dat;
