#!/bin/bash

# Get Pltlist
rm pltlist;
'/bin/ls' -1d plt????? >> pltlist;

# Job A 
jid1=$(sbatch --parsable DoProcess3DPartA)

# Job B
jid2=$(sbatch --dependency=afterok:${jid1} --parsable DoProcess3DPartBqslim)

# Job C
jid3=$(sbatch --dependency=afterok:${jid2} --parsable DoProcess3DPartCstream)

# Job D
sbatch --dependency=afterok:${jid3} DoProcess3DPartDstreamTubeStats 
