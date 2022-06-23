import numpy as np
import os
import matplotlib.pyplot as plt

pltlist = np.genfromtxt('pltlist',dtype='str')
nFiles = len(pltlist)
SC = ['FF','LP','LE','TE','TP','SP']
comps = ['H2_MassFraction','H2_ConsumptionRate','Temp'] #add phi here when done

pathLength = os.environ['NRK']
array = np.zeros((len(SC),len(comps),pathLength))
for plot in pltlist:
    filename=plot+'_prog_H2_'+os.environ['i_H2']+'_'+os.environ['NPTS']+'_streamBin'
    for i in range(len(SC)):
        for j in range(len(comps)):
            s_filename = filename+'/'+SC[i]+'_'+comps[j]+'.dat'
            localstream = np.loadtxt(s_filename)
            print(localstream)
            array[i,j] += localstream/nFiles
for i in range(len(SC)):
    for j in range(len(comps)):
        np.savetxt('PKZstreams/'+SC[i]+'_'+comps[j]+'.dat',array[i,j])
