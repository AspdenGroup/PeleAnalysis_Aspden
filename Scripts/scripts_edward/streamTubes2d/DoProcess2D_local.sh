#!/bin/bash
set -x

export ntasks=8
export isoVal=0.9
export fuel=H2

export fmf_loc=9 # H2 newPeleLM
export rr_loc=22 # H2 newPeleLM
export hr_loc=23 # H2 newPeleLM
#export fmf_loc=11 # CH4 newPeleLM
#export rr_loc=25 # CH4 newPeleLM
#export hr_loc=26 # CH4 newPeleLM
#export fmf_loc=60 # C12H26 newPeleLM
#export rr_loc=69 # C12H26 newPeleLM
#export hr_loc=70 # C12H26 newPeleLM

export NCELLSPERLF=16
export NPTS=100000
#export JHI=`echo "12*${NCELLSPERLF}" | bc -l`
#export NRK=`echo "2*${JHI} + 1" | bc -l`
#export HRK="0.25"
#export GROW=`echo "3*${NCELLSPERLF} + 2" | bc -l`
export JHI="120"
export NRK=`echo "2*${JHI} + 1" | bc -l`
export HRK="0.3"
    
./time.sh;
for PLT in plt?????;
do
    echo "--- running plotProg ---"
    mpirun -n ${ntasks} plotProg2d.gnu.MPI.ex inputs.process infile= ${PLT} fuelName=${fuel}

    echo "--- running grad ---"
    mpirun -n ${ntasks} grad2d.gnu.MPI.ex inputs.process infile= ${PLT}

    echo "--- Combining Plots ---"
    mpirun -n ${ntasks} combinePlts2d.gnu.MPI.ex inputs.process infileL= ${PLT}_prog \
	   infileR= ${PLT} \
	   compsL= 0 1  \
	   compsR= 0 1 ${hr_loc} ${rr_loc} \
	   outfile= ${PLT}_prog_xy

    echo "--- running curvature ---"
    mpirun -n ${ntasks} curvature2d.gnu.MPI.ex inputs.curvature infile= ${PLT}_prog_xy

    echo "--- Combining Plots ---"
    mpirun -n ${ntasks} combinePlts2d.gnu.MPI.ex inputs.process infileL= ${PLT}_gt \
	   infileR= ${PLT}_prog_xy_K \
	   compsL= 0 3 \
	   compsR= 0 3 7 10  \
	   outfile= ${PLT}_combine

    echo "--- running isosurface ---"
    mpirun -n ${ntasks} isosurface2d.gnu.MPI.ex inputs.process \
	   comps= 2 \
	   isoCompName= prog_${fuel} \
	   isoVal= ${isoVal} \
	   infile= ${PLT}_combine

    echo "--- running streamlines ---"
    # may cause issues
    mpirun -n ${ntasks} stream2d.gnu.MPI.ex inputs.process \
	   progressName= prog_${fuel} \
	   plotfile=${PLT}_combine \
	   isoFile=${PLT}_combine_prog_${fuel}_${isoVal}.mef \
	   streamFile=${PLT}_stream \
	   hRK=${HRK} nRKsteps=${NRK}

    echo "--- running sample streamlines ---"
    mpirun -n ${ntasks} sampleStreamlines2d.gnu.MPI.ex inputs.process \
	   nCompsPerPass=1 \
	   nGrow=40 \
	   plotfile=${PLT}_combine \
	   pathFile=${PLT}_stream \
	   streamSampleFile=${PLT}_stream_sample \
	   comps= 0 1 2 3 4 5

    echo "--- running streamtubestats ---"
    ./streamTubeStats2d.gnu.ex \
	verbose=1 \
	infile= ${PLT}_stream_sample \
	intComps= 6 peakComp= 4 6 avgComps= 7 8 
    
    echo "--- converting mef to dat ---"
    mpirun -n ${ntasks} surfMEFtoDATbasic2d.gnu.MPI.ex infile= ${PLT}_stream_sample_volInt.mef
    mpirun -n ${ntasks} surfMEFtoDAT2d.gnu.MPI.ex infile= ${PLT}_stream_sample_volInt.mef  
done
