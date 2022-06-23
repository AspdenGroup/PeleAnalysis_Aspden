import numpy as np
import os
import matplotlib.pyplot as plt

pltlist = np.genfromtxt('pltlist',dtype='str')
nFiles = len(pltlist)
SC = ['FF','LP','LE','TE','TP','SP']
comps = ['distfromseed','H2_MassFraction','H2_ConsumptionRate','Temp'] #add phi here when done

pathLength = int(os.environ['NRK'])
array = np.zeros((len(SC),len(comps),pathLength))
zonels = np.zeros(len(SC))
for plot in pltlist:
    filename=plot+'_prog_H2_'+os.environ['i_H2']+'_'+os.environ['NPTS']+'_streamBin'
    for i in range(len(SC)):
        zonels[i] += np.loadtxt(filename+'/'+SC[i]+'_ls.dat')/nFiles
        for j in range(len(comps)):
            s_filename = filename+'/'+SC[i]+'_'+comps[j]+'.dat'
            localstream = np.loadtxt(s_filename)
            array[i,j] += localstream/nFiles
for i in range(len(SC)):
    np.savetxt('properties/'+SC[i]+'_ls',[zonels[i]])
    for j in range(len(comps)):
        np.savetxt('PKZstreams/'+SC[i]+'_'+comps[j]+'.dat',array[i,j])

#plt.plot(array[0,0]/zonels[0],array[0,1])
#plt.show()
