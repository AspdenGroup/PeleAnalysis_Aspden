import numpy as np
import math
import pandas as pd

Is_3D = 0

filename = "int_prog.dat"
data = np.loadtxt(filename)
t = data[:,0]
V = data[:,1]
if (Is_3D == 1):
    r = ((V / math.pi) * (3 / 4))**(1 / 3)
else :
    r = np.sqrt(V / math.pi)
pd.DataFrame(r,t).to_csv('radius.csv',header=None)
