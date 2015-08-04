#!/usr/bin/python
"""
A script that creates a lut for input to the bash script that will run
a python script per entry in the lut. This can be modified to add a
distributed selection of parameters to test against.
"""
import numpy as np

sol_zen = np.arange(20., 90., 10.) # degrees
sat_zen = np.arange(0., 80., 10.)
rel_azi = np.arange(0., 195., 15.) # relative azimuth in GOME-2 geometry
alt = np.linspace(0.0, 5.0, 4, endpoint=True) # km AMSL
aot = np.linspace(0.0, 1.0, 4, endpoint=True)
wvc = np.concatenate((np.arange(0., 12., 6.), np.array([12.]),\
    np.arange(16., 28., 6.), np.arange(28., 36., 4.), np.arange(36., 52., 8.),\
    np.array([52.]), np.linspace(56., 80., 3, endpoint=True))) # kg/m2
press = np.linspace(500., 1100., 4, endpoint=True) # hPa

tot = len(sol_zen)*len(sat_zen)*len(rel_azi)*len(alt)*len(aot)*len(wvc)*\
    len(press)

large = []

for i, sz in enumerate(sol_zen):
    for j, vz in enumerate(sat_zen):
        for k, ra in enumerate(rel_azi):
            for l, al in enumerate(alt):
                for m, ao in enumerate(aot):
                    for n, wv in enumerate(wvc):
                        for o, pr in enumerate(press):
                            large.append(np.array([sz, vz, ra, al, ao, wv, pr]))

large = np.array(large)

fn = 'input_lut.dat'
header = 'sol_zen sat_zen rel_azi alt AOT WVC press'
np.savetxt(fn, large, fmt='%.3f', delimiter=' ', header=header, newline='\n')
