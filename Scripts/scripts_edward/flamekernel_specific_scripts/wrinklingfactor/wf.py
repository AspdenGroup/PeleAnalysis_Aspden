import numpy as np
import math
import pandas as pd

Is_3D = 0

filename = "int_prog.dat"
data = np.loadtxt(filename)
V = data[:,1]
# Effective Area
if (Is_3D == 1):
    r = ((V / math.pi) * (3 / 4))**(1 / 3)
else :
    r = np.sqrt(V / math.pi)
if (Is_3D == 1):
    A_eff = 4 * math.pi * r^2      # Surface Area
else :
    A_eff = 2 * math.pi * r       # Surface Area
# Actual Area
filename = "area.dat"
data = np.loadtxt(filename)
t = data[:,0]
A = data[:,1]
# Wrinkling Factor
wf = A / A_eff

pd.DataFrame(wf,t).to_csv('wrinkling_factor.csv',header=False)
