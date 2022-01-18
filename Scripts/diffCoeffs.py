import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
# Simulation parameters
phi = 0.7
atm = 1
p = atm*ct.one_atm
Tin = 300
width = 1  # m
loglevel = 0  # amount of diagnostic output (0 to 8)
reactants = {'NC12H26':(0.21/18.5)*phi, 'O2':0.21, 'N2':0.79}  # premixed gas composition
#IdealGasMix object used to compute mixture properties, set to the state of the upstream fuel-air mixture
gas = ct.Solution('dodecane_wang+Ar.cti')
gas.TPX = Tin, p, reactants

fuel_name = 'NC12H26'
index = np.argwhere(np.array(gas.species_names) == fuel_name)[0,0]
print('Mixture averaged diffusion coeffients')
print(gas.mix_diff_coeffs)
print('Thermal diffusion coefficient')
print(gas.thermal_conductivity/(gas.density*gas.cp_mass))
lewis_numbers = gas.thermal_conductivity/(gas.density*gas.mix_diff_coeffs*gas.cp_mass)
print('Species Lewis numbers')
print(lewis_numbers)
print('Fuel Lewis number')
print(lewis_numbers[index])
"""
f = ct.FreeFlame(gas,width=width)
f.set_refine_criteria(ratio=3, slope=0.005, curve=0.005)
f.max_grid_points = 10000
# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)
""" 
