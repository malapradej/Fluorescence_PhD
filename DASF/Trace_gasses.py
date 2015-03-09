#!/usr/bin/python
"""
This is a script to extract H2O total column information from 
GOME-2 level 2 products.
"""

import numpy as np
import matplotlib.pylab as plt
from fnmatch import fnmatch
import h5py
import os

root = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/H2O/atmos.caf.dlr.de/gome2a/offline/2007'
pattern = '*.HDF5'
Fns = []

for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            Fns.append(os.path.join(path, name))

Lat = []
Lon = []
TC_H2O = []
TC_H2O_err = []

nfiles = len(Fns)

for i, fn in enumerate(Fns):
    print 'file %d or %d.' % (i, nfiles)
    f = h5py.File(fn, 'r')
    
    # select ordering of gasses
    gas = 'H2O'
    gas_nr = np.where(f[u'META_DATA'][u'MainSpecies'][:]==gas)[0][0]
    
    # total column
    # units of kg/m2 (libRadtran MM)
    # flag needs to be 0 for correct retrieval
    tc_H2O_flag = f[u'DETAILED_RESULTS'][u'QualityFlags'][:,gas_nr]  # column 5
    index_flag = np.where(tc_H2O_flag==0)[0]          # index of correct retrievals
    tc_H2O = f["/TOTAL_COLUMNS"][u'H2O'][:]
    index_fill = np.where(tc_H2O!=-1)[0]              # fill value of -1
    tc_H2O_err = f["/TOTAL_COLUMNS"][u'H2O_Error'][:]           # % error
    
    # following reading it is decided that O2 can best be set using 
    # surface pressure. see atmospheric report I wrote....
    
    # effective slant col
    # needs to be converted to vertical col by eq(3) in ATBD
    # no molecular Ring correction needs to be applied see p.34 in ATBD
    # see notes on 4/2/15
    # this requires AMF for clear and cloudy sky, ghost column and 
    # intensity weighted cloud fraction
    # units of mol/cm2 (libRadtran CM_2)
    #esc_O2 = f[u'DETAILED_RESULTS'][u'H2O'][u'ESC_O2'][:]   
    #esc_O4 = f[u'DETAILED_RESULTS'][u'H2O'][u'ESC_O4'][:] # O4 not in libRadtran
    #iw_cf = f[u'DETAILED_RESULTS'][u'IntensityWeightedCloudFraction'][:,gas_nr]
    #ghost_c = f[u'DETAILED_RESULTS'][u'GhostColumn'][:,gas_nr]
    #amf_total = f[u'DETAILED_RESULTS'][u'AMFTotal'][:,gas_nr]
    #amf_cloud = f[u'DETAILED_RESULTS'][u'AMFToCloudTop'][:,gas_nr]
    #vcd_O2 = (esc_O2 + iw_cf * ghost_c * amf_cloud) / amf_total
    #vcd_O2 = esc_O2 / amf_total 
    
    # centre of pixel lat and lon
    # units of decimal degree
    lat = f[u'GEOLOCATION'][u'LatitudeCentre'][:]
    lon = f[u'GEOLOCATION'][u'LongitudeCentre'][:]
    
    # view mode in binary
    # must be nadir (nominal) and descending part of orbit
    # which is 100000000 or 256
    view_mode = f[u'GEOLOCATION'][u'ViewMode'][:]
    index_nadir = np.where(view_mode==256)[0] # index of nadir pixels
    
    index_good = list(set(index_flag.flatten()) & set(index_nadir.flatten())\
            & set(index_fill.flatten()))
    
    lat = lat[index_good]
    lon = lon[index_good]
    tc_H2O = tc_H2O[index_good]
    tc_H2O_err = tc_H2O_err[index_good]
    
    Lat.extend(list(lat))
    Lon.extend(list(lon))
    TC_H2O.extend(list(tc_H2O))
    TC_H2O_err.extend(list(tc_H2O_err))
    
    f.close()


plt.hexbin(Lat, TC_H2O, bins='log', cmap=plt.cm.Blues)
plt.xlabel('Latitude(S-, N+)', fontsize='large')
plt.ylabel(r'TOTAL COLUMN $H_2O (kg/m^2)$', fontsize='large')
#plt.plot(Lat, tc_H2O_err)
plt.show()

weights = np.ones_like(TC_H2O)/len(TC_H2O)
plt.hist(TC_H2O, bins=20, normed=0, weights=weights)
plt.xlabel(r'TOTAL COLUMN $H_2O (kg/m^2)$', fontsize='large')
plt.ylabel('fraction', fontsize='large')
plt.show()