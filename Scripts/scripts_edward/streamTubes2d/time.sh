#!/bin/bash

rm time.csv
for i in plt?????
do  
    sed '28q;d' $i/Header >> time.csv  ## H2 newPeleLM
    ##sed '43q;d' $i/Header >> time.csv  ## CH4 newPeleLM
    ##sed '75q;d' $i/Header >> time.csv  ## C12H26 newPeleLM
done
