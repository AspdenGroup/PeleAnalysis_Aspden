import numpy as np
import math
import pandas as pd

use_cantera = 1
Is_3D = 0

# Setup - only if use cantera #
phi = 0.4
atm = 40
Tin = 700
Is_LiDryer = 1
Is_drm19 = 0
Is_dodecane_wang = 0
###############################

if use_cantera == 0:
    YH2 = 100
    rho = 100
else:
    import cantera as ct    
    # Run Cantera ##
    print("Running Cantera")
    p = atm*ct.one_atm
    width = 1
    loglevel = 0
    if (Is_LiDryer == 1):
        reactants = {'H2':0.42*phi, 'O2':0.21, 'N2':0.79}
        gas = ct.Solution('LiDryer.cti')
        fuel_loc = 0
    if (Is_drm19 == 1):
        A = 1 / (0.105 * phi)
        reactants = {'CH4':1, 'O2':0.21*A, 'N2':0.79*A}
        gas = ct.Solution('drm19.cti')
        fuel_loc = 10
    if (Is_dodecane_wang == 1):
        A = 1 / (0.011351351 * phi)
        reactants = {'NC12H26':1, 'O2':0.21*A, 'N2':0.79*A}
        gas = ct.Solution('dodecane.cti')
        fuel_loc = 51
    gas.TPX = Tin, p, reactants
    f = ct.FreeFlame(gas,width=width)
    f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)
    f.transport_model = 'Mix'
    f.max_grid_points=10000
    f.solve(loglevel=loglevel, auto=True)
    
    YH2 = f.Y[fuel_loc,0]  
    rho = f.density[0] 

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
