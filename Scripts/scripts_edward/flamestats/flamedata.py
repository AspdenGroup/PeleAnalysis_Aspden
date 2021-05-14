import numpy as np
import math
import pandas as pd

# Setup #
use_cantera = 0
every100 = 0 # plot interval 0 = 10, 1 = 100

# setup # 
phi = 0.4
atm = 40
Tin = 700
Is_3D = 0
Is_LiDryer = 1
Is_drm19 = 0
Is_dodecane_wang = 0
#################################

# setup - only if not using cantera # 
if use_cantera == 0:
    Y_fuel = 0.01160255868152796       # H2 0.4phi-40atm-700K
    rho_fuel = 17.402904910349093      # H2 0.4phi-40atm-700K
    Tend = 1767.4265061553956          # H2 0.4phi-40atm-700K
    
#    Y_fuel = 0.055166072814902947     # CH4 1phi-40atm-700K
#    rho_fuel = 19.243251010010859     # CH4 1phi-40atm-700K
#    Tend = 2534.492660151539          # CH4 1phi-40atm-700K

#    Y_fuel = 0.062810985972055658     # C12H26 1phi-40atm-700K
#    rho_fuel = 21.196636610858601     # C12H26 1phi-40atm-700K
#    Tend = 2599.1153845893355         # C12H26 1phi-40atm-700K
#####################################

if use_cantera == 1:
    import cantera as ct
    # run cantera
    p = atm*ct.one_atm
    width = 1
    loglevel = 0
    if (Is_LiDryer == 1):
        reactants = {'H2':0.42*phi, 'O2':0.21, 'N2':0.79}
        gas = ct.Solution('~/src/PeleProduction_ISC2022/Submodules/PelePhysics/Support/Fuego/Mechanism/Models/LiDryer/LiDryer.cti')
        fuel_loc = 0
    if (Is_drm19 == 1):
        A = 1 / (0.105 * phi)
        reactants = {'CH4':1, 'O2':0.21*A, 'N2':0.79*A}
        gas = ct.Solution('~/src/PeleProduction_ISC2022/Submodules/PelePhysics/Support/Fuego/Mechanism/Models/drm19/drm19.cti')
        fuel_loc = 10
    if (Is_dodecane_wang == 1):
        A = 1 / (0.011351351 * phi)
        reactants = {'NC12H26':1, 'O2':0.21*A, 'N2':0.79*A}
        gas = ct.Solution('~/src/PeleProduction_ISC2022/Submodules/PelePhysics/Support/Fuego/Mechanism/Models/dodecane_wang/dodecane_reduced_mech_jetsurf.cti')
        fuel_loc = 51
    gas.TPX = Tin, p, reactants
    f = ct.FreeFlame(gas,width=width)
    f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)
    f.transport_model = 'Mix'
    f.max_grid_points=10000
    f.solve(loglevel=loglevel, auto=True)
    n_intervals = 100
    Y_fuel = f.Y[fuel_loc,0]
    rho_fuel = f.density[0]
    Tend = f.T[-1]
    
    print("Y_fuel = ",Y_fuel)
    print("rho_fuel = ",rho_fuel)
    print("f.T[-1] = ",f.T[-1])
    
   
# loads data and arrays
n_intervals = 100
time_data = np.loadtxt('time.csv')
n_pltfiles = len(time_data) -1 
plot = 0

mean_sloc = np.empty(n_pltfiles)
mean_lloc = np.empty(n_pltfiles)

sd_sloc = np.empty(n_pltfiles)
sd_lloc = np.empty(n_pltfiles)

t = np.empty(n_pltfiles)

for i in range(n_pltfiles):
    # sets up some more variables
    A_sloc_array = []
    sloc_array = []
    A_gt_array = [] 
    lloc_array = []
    
    # sets plotfile name
    if (every100 == 1):
        if (i < 1): 
            filename = 'plt00000'
        elif (i < 10):
            filename = 'plt00'+str(i*100)
        elif (i >= 10 and i < 100):
            filename = 'plt0'+str(i*100)
        else:
            filename = 'plt'+str(i*100)
        print(filename)
    else:
        if (i < 1): 
            filename = 'plt00000'
        elif (i < 10):
            filename = 'plt000'+str(i*10)
        elif (i >= 10 and i < 100):
            filename = 'plt00'+str(i*10)
        else:
            filename = 'plt0'+str(i*10)
        print(filename)       

    # Loads data files
    stream_data = np.loadtxt(filename+'_stream_sample_volInt_basic.dat')
    t[i] = time_data[i]
    for j in range(len(stream_data[:,1]))[::-1]:
        if (stream_data[j,11] != 1 or math.isnan(stream_data[j,3]) is True or math.isnan(stream_data[j,6]) is True or math.isnan( stream_data[j,7]) is True):
            stream_data = np.delete(stream_data, (j), axis=0)

    # ------------ Flame Speed ------------ #

    # finds area and the integral fuel consumption rate
    A = stream_data[:,3]
    fcr_int = stream_data[:,6]

    # calculating mean flame speed
    A_sloc_total = 0
    for j in range(len(fcr_int)):
        sloc = (1 / (rho_fuel * Y_fuel)) * fcr_int[j]
        sloc_array.append(sloc)
        A_sloc_array.append(A[j])
        A_sloc_total += A[j]
    mean_sloc[i] = np.dot(sloc_array,A_sloc_array) / A_sloc_total

    # calculating standard deviation
    xmu2 = 0
    for j in range(len(fcr_int)):
        sloc = (1 / (rho_fuel * Y_fuel)) * fcr_int[j]
        xmu2 = xmu2 + (sloc - mean_sloc[i])**2
    SD = (xmu2 / len(fcr_int))**0.5
    sd_sloc[i] = SD

    # ------------ Thermal Thickness ------------ #
    
    # calculating mean local thermal thickness
    A_gt_total = 0
    A_gt = stream_data[:,3]
    max_gt = stream_data[:,9]
    for j in range(len(max_gt)):
        lloc = (Tend - Tin) / max_gt[j]
        lloc_array.append(lloc)
        A_gt_array.append(A_gt[j])
        A_gt_total += A_gt[j]
    mean_lloc[i] = np.dot(lloc_array, A_gt_array) / A_gt_total

    # calculating standard deviation
    xmu2_lloc = 0
    for j in range(len(max_gt)):
        lloc = (Tend - Tin) / max_gt[j]
        xmu2_lloc = xmu2_lloc + (lloc - mean_lloc[i])**2
    SD_lloc = (xmu2_lloc / len(max_gt))**0.5
    sd_lloc[i] = SD_lloc

# ------------ Plot Data ------------ # 
    pd.DataFrame(mean_lloc[:],t[:]).to_csv('lloc_mean.csv',header=None)
    pd.DataFrame(sd_lloc,t).to_csv('lloc_mean_SD.csv',header=None)

    pd.DataFrame(mean_sloc[:],t[:]).to_csv('sloc_mean.csv',header=None)
    pd.DataFrame(sd_sloc,t).to_csv('sloc_mean_SD.csv',header=None)
