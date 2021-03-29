#!/bin/bash

set -x
export FINPLT=`/bin/ls -d1rt plt????? | tail -n 1`
#rm -r oldprocessing
#mv processing oldprocessing
#mkdir processing
#cd processing
#ln -s ~/src/PeleAnalysis/Src/*2d.gnu.ex .
#ln -s ~/src/PeleAnalysis_Curvature/PeleAnalysis/Src/curvature2d.gnu.ex .
#ln -s ~/src/PeleAnalysis/Src/ModelSpecificAnalysis/plotProg2d.gnu.ex .
#ln -s ../inputs inputs.process
#cp ../../inputs.curvature ../../PDFs.py ../../DoProcess2D .
#cd ..
for PLT in plt?????; do
    cd processing
    #mkdir ${PLT}_processing_data
    #cp DoProcess2D ${PLT}_processing_data
    #cp ../../PDFs.py .
    cd ${PLT}_processing_data
    #ln -s ../../${PLT} .
    ln -s ../../plt00000 .
    #ln -s ../*.ex .
    #ln -s ../inputs.* .
    #bash DoProcess2D ${PLT} 16
    if [ "${PLT}" == "${FINPLT}" ]; then
       ln -s ~/src/forceData/AmrDerive/AmrDeriveMakePPM2d.Linux.g++.gfortran.ex .
       ./AmrDeriveMakePPM2d.Linux.g++.gfortran.ex infile= ${PLT}_prog vars= prog_temp prog_H2 useminmax1= 0 1 useminmax2= 0 1
       ./AmrDeriveMakePPM2d.Linux.g++.gfortran.ex infile= ${PLT}_comb vars= H2_ConsumptionRate
       convert ${PLT}_prog/prog_temp.ppm temp.png
       convert ${PLT}_comb/H2_ConsumptionRate.ppm H2CR.png
       convert ${PLT}_prog/prog_H2.ppm H2.png
       ln -s ~/src/PeleAnalysis/Src/streamConvert2d.gnu.ex .
       ./streamConvert2d.gnu.ex inputs.process nComp= 1 infile= ${PLT}_stream_sample_H2_9
       cp ../../../plotPathsPic.py .
       python plotPathsPic.py
       convert -trim paths.png paths_trim.png
       cp paths_trim.png ../paths.png
    fi
       
    cd ../..
done



    
    
