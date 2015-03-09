#!/usr/bin/python
'''
Script to create Lookup tables using libRadtran of atmospheric path reflectance, 
double transmittance and spherical albedo. This script only creates the LUT 
which can be used to create the Emulator. 
Run the script with the parameters as arguments as follows:

'''
#Section below is only to import modules required and load the pickled dictionary file created in the DASF_GOME-2.py script. 

import os
import sys
import numpy as np
import pickle
import matplotlib.pylab as plt
import re
import subprocess
from scipy.ndimage.filters import gaussian_filter1d as gf
from scipy import interpolate
import collections
from multiprocessing import Process, Queue
import pdb

fn = 'GOME_xxx_1B_M02_20070125124939Z_20070125143057Z_R_O_20120209112017Z.p'
data_path = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/DASF/'
fn = os.path.join(data_path, fn)
fl = open(fn, 'rb')
ds = pickle.load(fl)

# create mask for part of globe we need data for

mask = np.logical_and(np.logical_and(np.array(ds['Lon']) < -52, np.array(ds['Lon']) > -61), \
                  np.logical_and(np.array(ds['Lat']) < 2, np.array(ds['Lat']) > -8))

# various functions used later on

def rel_azimuth(Sol_azi, Sat_azi):
    '''A function that returns the relative azimuth angle difference 
    between the sun and the satellite, where sun is at zero. This is
    relative to the libRadtran geometry. See notes on 4/3/15....
    '''
    Sol_azi = np.array(Sol_azi) + 180.
    Sol_azi = np.where(Sol_azi >= 360., Sol_azi - 360., Sol_azi)
    Sat_azi = np.array(Sat_azi)
    rel = Sat_azi - Sol_azi
    rel = np.where(rel < 0., rel + 360., rel)
    
    return rel

Lon = np.array(ds['Lon'])[mask]
Lat = np.array(ds['Lat'])[mask]
Sat_zen = np.array(ds['Sat_zen'])[mask]
Sat_azi = np.array(ds['Sat_azi'])[mask]
Sol_zen = np.array(ds['Sol_zen'])[mask]
Sol_azi = np.array(ds['Sol_azi'])[mask]
Rel_azi = rel_azimuth(Sol_azi, Sat_azi)
Ref_toa = np.array(ds['Ref_toa'])[mask]
Lam = np.array(ds['Spec'])[mask]
Alt = np.array(ds['Alt'])[mask]
CF = np.array(ds['Cloud_frac'])[mask]

# select only a small subset
indices = [int(i) for i in sys.argv[1:]]
#indices = [0, 3, 4, 8, 9, 13]
Sat_zen = [Sat_zen[i] for i in indices]
Sat_azi = [Sat_azi[i] for i in indices]
Sol_zen = [Sol_zen[i] for i in indices]
Sol_azi = [Sol_azi[i] for i in indices]
Rel_azi = [Rel_azi[i] for i in indices]
Alt = [Alt[i]/1000. for i in indices] # convert to km AMSL

#sol_zena = np.array([40.]) # degrees
#sat_zena = np.array([0., 10., 20., 30.])
#rel_azia = np.array([170., 340., 350.]) # relative azimuth in libRadtran geometry
#alta = np.array([0., 0.250]) # km AMSL
atma = ['tropics'] #np.arange(['tropics', 'midlatitude_summer', 'midlatitude_winter', 'subarctic_summer', 'subarctic_winter'])
aota = np.array([0., 0.2, 0.4, 0.6, 0.8, 1.0]) 
wvca = np.array([0., 10., 20., 30., 40., 50., 60., 70., 80.]) # kg/m2
pressa = np.array([900., 9050., 1000, 1050., 1100.]) # hPa
file_leaf = 'leaf_spectrum.txt'
wl_min, wl_max = (540., 760.)
w = np.genfromtxt(file_leaf)
wl = w[:,0]
w = np.sum(w[:,1:], axis=1)
index = np.logical_and(wl <= wl_max, wl >= wl_min)
wl = wl[index]
w = w[index]
wc = np.vstack((wl, w)).T
inp_template = 'data_files_path /usr/local/share/libRadtran/data/\
    \noutput_quantity transmittance \nmol_abs_param reptran \npseudospherical \
    \nrte_solver disort \naerosol_default \naerosol_species_file continental_average \
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran \
    \natmosphere_file %s\
    \naerosol_set_tau_at_wvl 550 %.6f \nmol_modify H2O %.6f MM \
    \npressure %.6f \nwavelength %.8f %.8f \
    \nsza %.6f \nphi0 0.0 \numu %.6f \nphi %.6f \
    \nalbedo %.6f \naltitude %.3f \nzout %s \noutput_user lambda %s'

# the functions used for each parameter such as spherical albedo etc.

def trans_double(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press, que):
    '''A function that calculates the 2-way upward and downward transmittance.
    '''
    t1_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, np.cos(sat_zen*np.pi/180.), rel_azi, 0.0, \
        alt, 'sur', 'eglo')
    t2_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sat_zen, np.cos(sat_zen*np.pi/180.), rel_azi, 0.0, \
        alt, 'sur', 'eglo')
    # the downward transmittance
    process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
        subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    t1_arr, err1 = process.communicate(input=t1_inp)
    t1_arr = re.split('[\n\s]+', t1_arr)[1:-1]
    t1_arr = np.array(map(float, t1_arr))
    t1_arr = np.reshape(t1_arr, (-1,2))
    lam = t1_arr[:,0]
    t1_arr = t1_arr[:,1]
    # take in consideration the cosine of solar zenith see eq(6.6) in manual
    t1_arr = t1_arr / np.cos(sol_zen*np.pi/180.)
    
    # the upward transmittance
    process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
        subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    t2_arr, err2 = process.communicate(input=t2_inp)
    t2_arr = re.split('[\n\s]+', t2_arr)[1:-1]
    t2_arr = np.array(map(float, t2_arr))
    t2_arr = np.reshape(t2_arr, (-1,2))
    lam = t2_arr[:,0]
    t2_arr = t2_arr[:,1]
    # take in consideration the cosine of solar zenith see eq(6.6) in manual
    t2_arr = t2_arr / np.cos(sat_zen*np.pi/180.)
    
    tt_arr = t1_arr * t2_arr
    
    que.put((lam, tt_arr))

def spher_alb(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press, que):
    '''A function that calculates the spherical albedo of the sky.
    '''
    sph_alb_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, np.cos(sat_zen*np.pi/180.), rel_azi, 1.0, \
        alt, 'sur', 'spher_alb')
    sph_alb_inp = sph_alb_inp + '\ndisort_spherical_albedo' # bug needs this at end...
    process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    sph_alb_arr, err4 = process.communicate(input=sph_alb_inp)
    sph_alb_arr = re.split('[\n\s]+', sph_alb_arr)[1:-1]
    sph_alb_arr = np.array(map(float, sph_alb_arr))
    sph_alb_arr = np.reshape(sph_alb_arr, (-1,2))
    lam = sph_alb_arr[:,0]
    sph_alb_arr = sph_alb_arr[:,1]
    
    que.put((lam, sph_alb_arr))

def atm_path_refl(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press, que):
    '''A function that calculates the atmospheric path reflectance also called
    intrinsic atmospheric reflectance.
    '''
    ref_atm_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, np.cos(sat_zen*np.pi/180.), rel_azi, 0.0, \
        alt, 'toa', 'uu')
    process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    ref_atm_arr, err3 = process.communicate(input=ref_atm_inp)
    ref_atm_arr = re.split('[\n\s]+', ref_atm_arr)[1:-1]
    ref_atm_arr = np.array(map(float, ref_atm_arr))
    ref_atm_arr = np.reshape(ref_atm_arr, (-1,2))
    lam = ref_atm_arr[:,0]
    ref_atm_arr = ref_atm_arr[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)
    
    que.put((lam, ref_atm_arr))

def app_refl(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press, surf_alb):
    '''A function that calculates the apparent reflectance at the sensor.
    '''
    lam, ref_atm_arr = atm_path_refl(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press)
    # interpolate surface albedo to same wavelength as libRadtran
    if isinstance(surf_alb, collections.Iterable):
        f = interpolate.interp1d(surf_alb[:,0], surf_alb[:,1])
        surf_alb = f(lam)
    tt_arr = trans_double(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press)[1]
    sph_alb_arr = spher_alb(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press)[1]
    app_refl = ref_atm_arr + tt_arr*surf_alb / (1 - sph_alb_arr*surf_alb)
    
    return (lam, app_refl)

def convol_spec(lam, spec, fwhm):
    '''Convolves a spectrum with a Gaussina ILS with conversion from a FWHM.
    '''
    inter = (lam[1] - lam[0] + lam[-1] - lam[-2]) / 2.
    std = fwhm / 2. / np.sqrt(2. * np.log(2.)) 
    pix_std = std / inter
    G = gf(spec, pix_std)
    
    return G

def LUT_select(sol_zena, sat_zena, rel_azia, alta, wl_min, \
        wl_max, atma, aota, wvca, pressa):
    '''A function that creates a look-up table for all the different iterable
    parameters provided. Version for only a select number of points.
    '''
    tots = len(sol_zena)*len(atma)*len(aota)*len(wvca)*len(pressa)
    ite = 0
    descr = ['lam', 'atm_path', 'dbl_trans', 'spher_alb', 'sol_zen', 'sat_zen', 'rel_azi',\
             'alt', 'atm', 'AOT', 'WVC', 'press']
    lut = {k: list() for k in descr}
    for sz, vz, ra, al, in zip(sol_zena, sat_zena, rel_azia, alta):
        for at in atma:
            for ao in aota:
                for wv in wvca:
                    for pr in pressa:
                        ite += 1
                        sys.stdout.write('\r%d/%d' % (ite, tots))
                        sys.stdout.flush()
                        que1 = Queue()
                        que2 = Queue()
                        que3 = Queue()
                        p1 = Process(target=atm_path_refl, args=(sz, vz, ra,  al, wl_min, \
                            wl_max, at, ao, wv, pr, que1))
                        p1.start()
                        p2 = Process(target=trans_double, args=(sz, vz, ra,  al, wl_min, \
                            wl_max, at, ao, wv, pr, que2))
                        p2.start()
                        p3 = Process(target=spher_alb, args=(sz, vz, ra,  al, wl_min, \
                            wl_max, at, ao, wv, pr, que3))
                        p3.start()
                        lam, apr = que1.get()
                        p1.join()
                        lam, dtt = que2.get()
                        p2.join()
                        lam, sal = que3.get()
                        p3.join()
                        vals = [lam, apr, dtt, sal, sz, vz, ra, al, at, ao, wv, pr]
                        for k, v in zip(descr, vals):
                            lut[k].append(v)
    sys.stdout.write("\n")
    return lut

# <codecell>

lut = LUT_select(Sol_zen, Sat_zen, Rel_azi, Alt, wl_min, \
        wl_max, atma, aota, wvca, pressa)

# the LUT fn below is commented out and will be the template for the final one....

'''def LUT(sol_zen, sat_zen, rel_azi, alta, wl_min, \
        wl_max, atma, aota, wvca, pressa):
    A function that creates a look-up table for all the different iterable
    parameters provided.
    tots = len(sol_zena)*len(sat_zena)*len(rel_azia)\
        *len(alta)*len(atma)*len(aota)*len(wvca)*len(pressa)
    ite = 0
    descr = ['lam', 'atm_path', 'dbl_trans', 'spher_alb', 'sol_zen', 'sat_zen', 'rel_azi',\
             'alt', 'atm', 'AOT', 'WVC', 'press']
    lut = {k: list() for k in descr}
    for sz in sol_zena:
        for vz in sat_zena:
            for ra in rel_azia:
                for al in alta:
                    for at in atma:
                        for ao in aota:
                            for wv in wvca:
                                for pr in pressa:
                                    ite += 1
                                    sys.stdout.write('\r%d/%d' % (ite, tots))
                                    sys.stdout.flush()
                                    que1 = Queue()
                                    que2 = Queue()
                                    que3 = Queue()
                                    p1 = Process(target=atm_path_refl, args=(sz, vz, ra,  al, wl_min, \
    wl_max, at, ao, wv, pr, que1))
                                    p1.start()
                                    p2 = Process(target=trans_double, args=(sz, vz, ra,  al, wl_min, \
    wl_max, at, ao, wv, pr, que2))
                                    p2.start()
                                    p3 = Process(target=spher_alb, args=(sz, vz, ra,  al, wl_min, \
    wl_max, at, ao, wv, pr, que3))
                                    p3.start()
                                    lam, apr = que1.get()
                                    p1.join()
                                    lam, dtt = que2.get()
                                    p2.join()
                                    lam, sal = que3.get()
                                    p3.join()
                                    vals = [lam, apr, dtt, sal, sz, vz, ra, al, at, ao, wv, pr]
                                    for k, v in zip(descr, vals):
                                        lut[k].append(v)
    sys.stdout.write("\n")
    return lut'''

'''
lut = LUT(sol_zena, sat_zena, rel_azia, alta, wl_min, \
        wl_max, atma, aota, wvca, pressa)
'''
# save the data to disk...
exten = ''
for i in sys.argv[1:]:
    exten += '_%s' % i
fn = 'LUT%s.p' % exten
fl = open(fn, 'wb')
pickle.dump(lut, fl)
fl.close()

# 
'''
fn = 'LUT.p'
fl = open(fn, 'rb')
lut_test = pickle.load(fl)
fl.close()
'''
# 
'''
lut_test.keys()
'''
# 


