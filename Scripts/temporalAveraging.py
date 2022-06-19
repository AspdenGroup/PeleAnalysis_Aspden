import numpy as np
import os
import matplotlib.pyplot as plt

pltlist = np.genfromtxt('pltlist',dtype='str')
ss = 0
ls = 0
Ebar = 0
nFiles = len(pltlist)
#SC = ['FF','LP','LE','TE','TP','SP']
#comps = ['H2_MassFraction','H2_ConsumptionRate','Temp']

#pathLength = os.environ['NRK']
#array = np.zeros((len(SC),len(comps),pathLength))
for plot in pltlist:
    filename=plot+'_prog_H2_'+os.environ['i_H2']+'_'+os.environ['NPTS']+'_streamBin'
    characteristics = np.loadtxt(filename+'/characteristics.dat')
    ls += characteristics[0]/nFiles; 
    ss += characteristics[1]/nFiles;
    Ebar += characteristics[2]/nFiles;
#    for i in range(len(SC)):
#        for j in range(len(comps)):
#            s_filename = filename+'/'+SC[i]+'_'+comps[j]+'.dat'
#            localstream = np.loadtxt(s_filename)
#            print(localstream)
#            array[i,j] += localstream/nFiles
            
#print(ss)
#print(ls)
#print(Ebar)
np.savetxt('properties/ss',[ss])
np.savetxt('properties/ls',[ls])
np.savetxt('properties/Ebar',[Ebar])

#for j in range(len(comps)):
#plt.figure(1)
#for i in range(len(SC)):
#    plt.plot(array[i,2],array[i,0],label=SC[i])
#plt.legend()
#plt.show()
