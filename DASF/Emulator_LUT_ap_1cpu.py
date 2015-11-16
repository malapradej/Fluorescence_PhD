#!~/venv/bin/python
'''
LUT entry for libRadtran toa relfectance. This notebook is for the final LUT
which will be used to create the fast forward model used in surface reflectance
retrievals. A random selection of distributed atmospheric parameters will be
used as inputs to the Emulator.
The distribution parameters are included in the script. To run the script
simply run ./Emulator_LUT.py <sol_zen> <sat_zen> <rel_azi> <alt> <aot> <wvc>
    <press> <a> <p> <filename>
The script writes the LUT entry to a pickled dictionary with filename:
lut_<sol_zen>_<sat_zen>_<rel_azi>_<alt>_<aot>_<wvc>_<press>_<a>_<p>.p
rel_azi is in GOME-2 converntion where the solar azimuth is from evalution
location towards the sun, not in libRadtran which is from the direction of
the sun to the evalution location.
'''

import sys
import os
import errno
import numpy as np
import cPickle
import re
import subprocess
from scipy import signal
from scipy.interpolate import interp1d
import pdb

sol_zen = np.float(sys.argv[1]) #np.arange(20., 90., 10.) # degrees
sat_zen = np.float(sys.argv[2]) #np.arange(0., 80., 10.)
rel_azi = np.float(sys.argv[3]) #np.arange(0., 195., 15.) # relative azimuth in GOME-2 geometry
alt = np.float(sys.argv[4]) #np.array([0., 0.5, 2., 5.]) # km AMSL
atm = 'US-standard' #np.arange(['tropics', 'midlatitude_summer', 'midlatitude_winter', 'subarctic_summer', 'subarctic_winter'])
aot = np.float(sys.argv[5]) #np.linspace(0.0, 1.0, 4, endpoint=True)
wvc = np.float(sys.argv[6]) #np.concatenate((np.arange(0., 12., 6.), np.array([12.]),\
    #np.arange(16., 28., 6.), np.arange(28., 36., 4.), np.arange(36., 52., 8.),\
    #np.array([52.]), np.linspace(56., 80., 2, endpoint=True))) # kg/m2
press = np.float(sys.argv[7]) #np.linspace(500., 1100., 4, endpoint=True) # hPa
a = np.float(sys.argv[8])
p = np.float(sys.argv[9]) # recollision probability
wl_min, wl_max = (540., 800.)

inp_template = 'data_files_path /usr/local/share/libRadtran/data/\
    \noutput_quantity transmittance \nmol_abs_param reptran \npseudospherical \
    \nrte_solver disort \naerosol_default \naerosol_species_file continental_average \
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran \
    \natmosphere_file %s\
    \naerosol_set_tau_at_wvl 550 %.6f \nmol_modify H2O %.6f MM \
    \npressure %.6f \nwavelength %.8f %.8f \
    \nsza %.6f \nphi0 0.0 \numu %.6f \nphi %.6f \
    \nalbedo %.6f \naltitude %.3f \nzout %s \noutput_user lambda %s'

fwhm = 0.4931 # fwhm of slit function nm
w_slit = 2 # width of slit function in nm

# load the leaf albedo data
leaf = np.loadtxt('leaf_spectrum.txt').T
leafLambda = leaf[0]
leafAlbedo = leaf[1] + leaf[2]

# reduce the leaf lam and albedo data size (speeds up interpolation...)
xleafFull = leafLambda.copy()
leafLambda = xleafFull[(xleafFull>=wl_min) * (xleafFull<=wl_max)]
leafAlbedo = leafAlbedo[(xleafFull>=wl_min) * (xleafFull<=wl_max)]

def slit_smooth(lam, spectrum, fwhm, width):
    '''
    Function that applies the gaussian slit function with known FWHM to spectrum.
    It will replace the boundaries of the smoothed spectrum which are affected by
    zero padding with the original spectrum.
    Input:
    lam - wavelenght array in nm
    spectrum - the spectrum to convolve
    fwhm - the full width half maximum of the gaussian in nm
    width - the width of the slit function in nm
    Output:
    convol - smoothed spectrum
    '''
    resol = np.average(np.diff(lam)) # average spectral resol of model in nm
    intervals = np.ceil((width / resol)+1) # the number of intervals in the width of slit
    slit = signal.gaussian(intervals, fwhm/2.3548201*intervals/width)
    convol = signal.fftconvolve(spectrum, slit/np.sum(slit), mode='same')
    hwidth = (intervals-1)/2
    convol[:intervals] = spectrum[:intervals]
    convol[-intervals-1:] = spectrum[-intervals-1:]
    return convol

def rel_azimuth(rel_azi):
    '''A function that returns the relative azimuth angle difference
    between the sun and the satellite, where sun is at zero. This is
    relative to the libRadtran geometry. See notes on 4/3/15....
    '''
    # GOME-2 convention
    Sol_azi = 0.
    Sat_azi = rel_azi
    # libRadtran convention
    Sol_azi = Sol_azi + 180.
    if Sol_azi >= 360.:
        Sol_azi -= 360.
    rel = Sat_azi - Sol_azi
    if rel < 0.:
        rel += 360.
#    Sol_azi = np.array(Sol_azi) + 180.
#    Sol_azi = np.where(Sol_azi >= 360., Sol_azi - 360., Sol_azi)
#    Sat_azi = np.array(Sat_azi)
#    rel = Sat_azi - Sol_azi
#    rel = np.where(rel < 0., rel + 360., rel)
    return rel

# the functions used for each parameter such as spherical albedo etc.

def trans_double(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
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

    return (lam, tt_arr)

def spher_alb(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
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

    return (lam, sph_alb_arr)

def atm_path_refl(sol_zen, sat_zen, rel_azi,  alt, wl_min, \
        wl_max, atm, aot, wvc, press):
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

    return (lam, ref_atm_arr)

def surf_refl(a, p, leafW):
    '''A function that calculated the surface reflectance at the top of canopy
    based on the p-theory BRF calculation using
    a = canopy interceptance * directional gap density
    and
    BRF = a * leafW / (1 - p * leafW)
    where p is recollision probability and leafW the single scattering albedo of
    the leaf.
    '''
    refl = a * leafW / (1. - p * leafW)
    return refl

def toa_refl(apr, dtt, sal, srefl):
    '''A function that returns the toa reflectance based on surface reflectance
    and atmospheric variables.
    '''
    refl = apr + dtt * srefl / (1 - sal * srefl)
    return refl

def LUT_select(sz, vz, ra, al, wl_min, wl_max, at, ao, wv, pr, a, p):
    '''A function that creates a look-up table entry based on the parameters
    provided.
    '''
    descr = ['lam', 'atm_path', 'dbl_trans', 'spher_alb', 'surf_refl',\
        'toa_refl', 'sol_zen', 'sat_zen', 'rel_azi', 'alt', 'AOT', 'WVC',\
        'press', 'a', 'p']
    lam, apr = atm_path_refl(sz, vz, ra,  al, wl_min, \
        wl_max, at, ao, wv, pr)
    apr = slit_smooth(lam, apr, fwhm, w_slit)
    lam, dtt = trans_double(sz, vz, ra,  al, wl_min, \
        wl_max, at, ao, wv, pr)
    dtt = slit_smooth(lam, dtt, fwhm, w_slit)
    lam, sal = spher_alb(sz, vz, ra,  al, wl_min, \
        wl_max, at, ao, wv, pr)
    sal = slit_smooth(lam, sal, fwhm, w_slit)
    # leaf albedo interpolated to obs lam
    leafW = (interp1d(leafLambda, leafAlbedo, kind='linear'))(lam)
    srefl = surf_refl(a, p, leafW)
    trefl = toa_refl(apr, dtt, sal, srefl)
    vals = [lam, apr, dtt, sal, srefl, trefl, sz, vz, rel_azi, al, ao, wv, pr,\
        a, p]
    lut = dict(zip(descr, vals))
    return lut

lrt_azi = rel_azimuth(rel_azi)
lut = LUT_select(sol_zen, sat_zen, lrt_azi, alt, wl_min, \
        wl_max, atm, aot, wvc, press, a, p)

# save the data to disk...
root = '/home/malapradej/Documents/PhD_UCL/Data/LUT/test_ap_DASF'
path = os.path.join(root, sys.argv[1], sys.argv[2], sys.argv[3])
try:
    os.makedirs(path)
except OSError:
    if not os.path.isdir(path):
    	raise
fn = '%s_%s_%s_%s_%s_%s.p' % (sys.argv[4], sys.argv[5], sys.argv[6],\
    sys.argv[7], sys.argv[8], sys.argv[9])
fn = os.path.join(path, fn)
fl = open(fn, 'wb')
cPickle.dump(lut, fl, -1) # use highest protocol binary format
fl.close()
