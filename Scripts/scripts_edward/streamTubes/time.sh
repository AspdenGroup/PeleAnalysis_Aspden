#!/bin/bash

rm time.csv
for i in plt?????
do  
    sed '28q;d' $i/Header >> time.csv  ## Hydrogen
    ##sed '43q;d' $i/Header >> time.csv  ## Methane
    ##sed '78q;d' $i/Header >> time.csv  ## Dodecane
done
