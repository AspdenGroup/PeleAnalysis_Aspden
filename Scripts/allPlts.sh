#!/bin/bash

set -x
export FINPLT=`/bin/ls -d1rt plt????? | tail -n 1`

mkdir processing
cd processing
ln -s ~/src/PeleAnalysis_Aspden/Src/*2d.gnu.ex .
ln -s ../inputs inputs.process
cp ../../inputs.curvature ../../DoProcess2D .
cd ..
for PLT in plt?????; do
    cd processing
    mkdir ${PLT}_processing_data
    cp DoProcess2D ${PLT}_processing_data
    cd ${PLT}_processing_data
    ln -s ../../${PLT} .
    ln -s ../../plt00000 .
    ln -s ../*.ex .
    ln -s ../inputs.* .
    bash DoProcess2D ${PLT} 16   
    cd ../..
done



    
    
