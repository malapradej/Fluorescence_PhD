#!/usr/bin/python
"""
A script to calculate the combined upward and downward transmittance, 
spherical albedo and atmospheric path reflectance given the input values of
atmosphere, AOT, WVC, altitude, surface pressure, solar zenith, solar 
azimuth, view zenith, view azimuth, day-of-year....
"""

import subprocess
import matplotlib.pylab as plt
import numpy as np
import re
from scipy.ndimage.filters import gaussian_filter1d as gf
from scipy import interpolate
import collections
import pickle
import pdb

sol_zen = 60.
sol_azi = 90.
sat_zen = 40.
sat_azi = 270.
alt = 0.
wl_min = 400.
wl_max = 800.
atm = 'tropics'
aot = 0.12
wvc = 9.7
press = 1013.25
fwhm = 0.4931

surf_alb = 0.5
file_leaf = 'leaf_spectrum.txt'
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
    \nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
    \nalbedo %.6f \naltitude %.3f \nzout %s \noutput_user lambda %s'

def trans_double(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
    '''A function that calculates the 2-way upward and downward transmittance.
    '''
    t1_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.), sat_azi, 0.0, \
        alt, 'sur', 'eglo')
    t2_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sat_zen, sat_azi, np.cos(sat_zen*np.pi/180.), sat_azi, 0.0, \
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
    
    return (lam, tt_arr)

def spher_alb(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
    '''A function that calculates the spherical albedo of the sky.
    '''
    sph_alb_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.), sat_azi, 1.0, \
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
    
    return (lam, sph_alb_arr)

def atm_path_refl(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
    '''A function that calculates the atmospheric path reflectance also called
    intrinsic atmospheric reflectance.
    '''
    ref_atm_inp = inp_template % (atm, aot, wvc, press, wl_min,\
        wl_max, sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.), sat_azi, 0.0, \
        alt, 'toa', 'uu')
    process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    ref_atm_arr, err3 = process.communicate(input=ref_atm_inp)
    ref_atm_arr = re.split('[\n\s]+', ref_atm_arr)[1:-1]
    ref_atm_arr = np.array(map(float, ref_atm_arr))
    ref_atm_arr = np.reshape(ref_atm_arr, (-1,2))
    lam = ref_atm_arr[:,0]
    ref_atm_arr = ref_atm_arr[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)
    
    return (lam, ref_atm_arr)

def app_refl(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press, surf_alb):
    '''A function that calculates the apparent reflectance at the sensor.
    '''
    lam, ref_atm_arr = atm_path_refl(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press)
    # interpolate surface albedo to same wavelength as libRadtran
    if isinstance(surf_alb, collections.Iterable):
        f = interpolate.interp1d(surf_alb[:,0], surf_alb[:,1])
        surf_alb = f(lam)
    tt_arr = trans_double(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press)[1]
    sph_alb_arr = spher_alb(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
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
    

def voigt(fwhm, x0, x):
    '''A pseudo-voigt approximation of a voigt function for the ILS.
    '''
    eta = 0.5 # the weight of the lorentzian
    G = np.exp(-np.log(2.) * ((x - x0)/fwhm/2.)**2)
    L = 1/ (1. + ((x - x0)/fwhm/2.)**2)
    V = eta * L + (1. - eta) * G
    
    return V

def gauss(fwhm, x0, x):
    '''A gaussian function based on the fwhm for the ILS.
    '''
    std = fwhm / 2. / np.sqrt(2. * np.log(2.))
    G = 1. / std / np.sqrt(2.*np.pi) * np.exp( - 0.5*((x - x0)/std)**2.)
    integ = 1. / std / np.sqrt(2.*np.pi) * std * 2. * np.pi
    G = G/integ
    
    return G

def BRF(rhoi0, p, wc):
    '''A function that calculates the BRF given the directional gap density,
    recollision probability, and leaf ss albedo.
    '''
    w = wc[:,1]
    wl = wc[:,0]
    brf = w * rhoi0 / (1. - w * p)
    temp = np.vstack((wl, brf)).T
    return temp

def LUT(sol_zena, sol_azia, sat_zena, sat_azia, alta, wl_min, \
        wl_max, atma, aota, wvca, pressa, rhoi0a, pa):
    '''A function that creates a look-up table for all the different iterable
    parameters provided.
    '''
    tots = len(sol_zena)*len(sol_azia)*len(sat_zena)*len(sat_azia)\
        *len(alta)*len(atma)*len(aota)*len(wvca)*len(pressa)*len(rhoi0a)\
        *len(pa)
    ite = 0
    descr = ['lam', 'app_refl', 'BRF', 'sol_zen', 'sol_azi', 'sat_zen', 'sat_azi',\
        'alt', 'atm', 'AOT', 'WVC', 'press', 'rhoi0', 'p']
    lut = []
    for sz in sol_zena:
        for sa in sol_azia:
            for vz in sat_zena:
                for va in sat_azia:
                    for al in alta:
                        for at in atma:
                            for ao in aota:
                                for wv in wvca:
                                    for pr in pressa:
                                        for rh in rhoi0a:
                                            for pp in pa:
                                                ite += 1
                                                print '%d/%d\n' % (ite, tots)
                                                brf = BRF(rh, pp, wc)
                                                lam, app_refl_arr = \
    app_refl(sz, sa, vz, va, al, wl_min, wl_max, at, ao, wv, pr, brf)
                                                lut.append([lam, \
    app_refl_arr, brf, sz, sa, vz, va, al, at, ao, wv, pr, rh, pp])
    
    return (descr, lut)

sol_zena = np.array([0., 60.]) #np.linspace(0., 80./180.*np.pi, 9)
sol_azia = np.array([0., 90.]) #np.linspace(0., 2.*np.pi, 36, endpoint=False)
sat_zena = np.array([0., 60.]) #np.linspace(0., 80./180.*np.pi, 9)
sat_azia = np.array([0.]) #np.linspace(0., 2.*np.pi, 36, endpoint=False)
alta = np.array([0., 4.]) #np.array([0., 0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 6., 7.])
atma = ['tropics', 'subarctic_winter'] #np.arange(['tropics', 'midlatitude_summer', 'midlatitude_winter', 'subarctic_summer', 'subarctic_winter'])
aota = np.array([0., 0.2]) #np.array([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,\
    #2., 4., 6.])
wvca = np.array([0., 30.]) #np.array([0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 25., 30., \
    #40., 50., 60., 80.])
pressa = np.array([950., 1100.]) #np.arange(500., 1100., 50.)
rhoi0a = np.array([0., 0.5, 1.]) #np.linspace(0., np.pi, 10)
pa = np.array([0., 0.5]) #np.linspace(0., 1., 11)

lut = LUT(sol_zena, sol_azia, sat_zena, sat_azia, alta, wl_min, \
        wl_max, atma, aota, wvca, pressa, rhoi0a, pa)

fl = open('LUT.p', 'wb')
pickle.dump(lut, fl)
fl.close()

def plotcomb(lut):
    '''A function that plots the desired fields for specific records.
    '''
    fixed_fields = ['sol_azi', 'sol_zen', 'sat_zen', 'sat_azi',\
        'alt', 'press', 'p', 'rhoi0']
    fixed_values = [0.0, 0.0, 0.0, 0.0, 0.0, 1050, 0.66, 0.5]
    print_fields = ['AOT', 'WVC', ]
    print_index = []
    for i in print_fields:
        print_index.append(lut[0].index(i))
    xy_pair = ['BRF', 'app_refl']
    fixed_index = []
    for i in fixed_fields:
        fixed_index.append(lut[0].index(i))
    #var_index = []
    #for i in var_fields:
    #    var_index.append(lut[0].index(i))
    rec_index = []
    for i, j in zip(fixed_index, fixed_values):
        temp_index = []
        for k, l in enumerate(lut[1][:]):
            if l[i] == j:
                temp_index.append(k)
        rec_index.append(temp_index)
    rec_index = set.intersection(*(set(i) for i in rec_index))
    for i in xy_pair:
        field = lut[0].index(i)
        if i == 'BRF':
            xs = [lut[1][j][field][:,0] for j in rec_index]
            ys = [lut[1][j][field][:,1] for j in rec_index]
            for j, (x, y, r) in enumerate(zip(xs, ys, rec_index)):
                st = ''
                for k, l in zip(print_fields, print_index):
                    st += '%s %.2f ' %(k, lut[1][r][l])
                label = '%s %s' %(i, st)
                plt.plot(x, y, label=label)
        else:
            lam_field = lut[0].index('lam')
            xs = [lut[1][j][lam_field] for j in rec_index]
            ys = [lut[1][j][field] for j in rec_index]
            for j, (x, y, r) in enumerate(zip(xs, ys, rec_index)):
                st = ''
                for k, l in zip(print_fields, print_index):
                    st += '%s %.2f ' %(k, lut[1][r][l])
                label = '%s %s' %(i, st)
                plt.plot(x, y, label=label)
    plt.xlabel('wavelength (nm)', fontsize='medium')
    plt.ylabel('Reflectance', fontsize='medium')
    plt.legend(loc='upper left', fontsize='small')
    plt.show()



'''
lam, app_refl_arr = app_refl(sol_zen, sol_azi, sat_zen, sat_azi,  alt, wl_min, \
    wl_max, atm, aot, wvc, press, surf_alb)

plt.plot(lam, app_refl_arr, 'r-')
plt.show()
'''