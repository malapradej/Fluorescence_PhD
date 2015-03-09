#!/usr/bin/python
"""
A script to extract aerosol optical thickness from SeaWIFS Aerosol product.
"""
import os
import glob
import numpy as np
import scipy.misc as misc
import matplotlib.pylab as plt
import gdal

root = '/home/malapradej/Documents/PhD_UCL/Data/SeaWIFS/Aerosols'
pattern = '*.h5'
crit = os.path.join(root, pattern)
Fns = glob.glob(crit)

crit_AOT = 'HDF5:\"%s\"://aerosol_optical_thickness_land'
crit_AOT_std = 'HDF5:\"%s\"://aerosol_optical_thickness_stddev_land'
crit_AOT_max = 'HDF5:\"%s\"://diagnostics/aerosol_optical_thickness_maximum_land'
crit_AOT_min = 'HDF5:\"%s\"://diagnostics/aerosol_optical_thickness_minimum_land'

Lat = []
Lon = []
AOT = []
AOT_std = []
AOT_max = []
AOT_min = []

for fn in Fns:
    dsn_aot = crit_AOT % fn
    ds_AOT = gdal.Open(dsn_aot)
    arr_AOT = ds_AOT.ReadAsArray()[2]
    arr_AOT = np.where(arr_AOT==-999, np.nan, arr_AOT)
    AOT.append(arr_AOT)
    
    dsn_aot_std = crit_AOT_std % fn
    ds_AOT_std = gdal.Open(dsn_aot_std)
    arr_AOT_std = ds_AOT.ReadAsArray()[2]
    arr_AOT_std = np.where(arr_AOT_std==-999, np.nan, arr_AOT_std)
    AOT_std.append(arr_AOT_std)
    
    dsn_aot_max = crit_AOT_max % fn
    ds_AOT_max = gdal.Open(dsn_aot_max)
    arr_AOT_max = ds_AOT_max.ReadAsArray()[2]
    arr_AOT_max = np.where(arr_AOT_max==-999, np.nan, arr_AOT_max)
    AOT_max.append(arr_AOT_max)
    
    dsn_aot_min = crit_AOT_min % fn
    ds_AOT_min = gdal.Open(dsn_aot_min)
    arr_AOT_min = ds_AOT_min.ReadAsArray()[2]
    arr_AOT_min = np.where(arr_AOT_min==-999, np.nan, arr_AOT_min)
    AOT_min.append(arr_AOT_min)


AOT_max_all = np.nanmax(AOT_max, axis=0)

rows, cols = np.shape(AOT_max_all)
x = np.linspace(-180, 180, cols+1)
y = np.linspace(-90, 90, rows+1)
xx, yy = np.meshgrid(x, y)
masked_arr = np.ma.array(AOT_max_all, mask=np.isnan(AOT_max_all))
plt.pcolormesh(xx, yy, masked_arr)
#cmap = plt.cm.Reds
#plt.imshow(AOT_max_all, interpolation='none', vmin=0., vmax=7., cmap=cmap,\
#        origin='lower')

ds_out = gdal.GetDriverByName('GTiff').Create('./figures/AOT_max_map.tif', \
        cols, rows, 1, gdal.GDT_Float32)
band_out = ds_out.GetRasterBand(1)
band_out.WriteArray(masked_arr[::-1,:])
band_out.FlushCache()
ds_out.SetGeoTransform((-180.0, 0.5, 0.0, 90.0, 0.0, -0.5))
#ds_out.SetProjection(ds_AOT.GetProjection())
ds_out = None

lat_AOT_max = np.nanmax(AOT_max_all,axis=1)
lat_AOT_min = np.nanmin(np.nanmin(AOT_min, axis=0),axis=1)

lat_bins = np.arange(-90, 90, 0.5)
'''
plt.barh(lat_bins, lat_AOT_max)
plt.xlabel('maximum AOT')
plt.ylabel('Latitude (S-, N+)')
'''
#plt.show()