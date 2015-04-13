#!/usr/bin/python
"""
A script to calculate libRadtran transmittances, spherical albedo, intrinsic atmopsheric
reflectance and apparent reflectance to compare with 6S.
"""

import subprocess
import matplotlib.pylab as plt
import numpy as np
import re

plt.cla()

sol_zen = 40.75
sol_azi = 0.
sat_zen = 31.23
sat_azi = 342.24
alt = 0.15
wl_min = 540.
wl_max = 755.
resol = 2.
wvc = 40.
aot = 0.2
press = 950.

surf_refl = 0.1

inp_default = 'data_files_path /usr/local/share/libRadtran/data/\n\
    output_quantity transmittance \nmol_abs_param reptran\n\
    atmosphere_file US-standard \nrte_solver disort \naerosol_default \
    \naerosol_species_file continental_average \
    \naerosol_set_tau_at_wvl 550 %.3f \nmol_modify H2O %.3f MM \
    \nwavelength 400.0 1400.0 \npressure %.3f\
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran'\
    %(aot, wvc, press)

inp_t = '\nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
    \nalbedo %.1f \naltitude %.3f \nzout %s \nwavelength %.8f %.8f \
    \noutput_user lambda %s'
    
inp_t1 = inp_t % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'sur', wl_min, wl_max, 'eglo')
inp_t2 = inp_t % (sat_zen, sat_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'sur', wl_min, wl_max, 'eglo') #uu
inp_ref_atm = inp_t % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'toa', wl_min, wl_max, 'uu')
inp_sph_alb = inp_t % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 1.0, alt/1000., 'sur', wl_min, wl_max, 'spher_alb')
inp_sph_alb = inp_sph_alb + '\ndisort_spherical_albedo' # bug needs this at end...
inp_app_refl_true = inp_t % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, surf_refl, alt/1000., 'toa', wl_min, wl_max, 'uu')
inp_t1 = inp_default + inp_t1
inp_t2 = inp_default + inp_t2
inp_ref_atm = inp_default + inp_ref_atm
inp_sph_alb = inp_default + inp_sph_alb
inp_app_refl_true = inp_default + inp_app_refl_true

# these are the shell pipes needed to run the command line uvspec program
# the downward transmittance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
t1_arr, err1 = process.communicate(input=inp_t1)
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
t2_arr, err2 = process.communicate(input=inp_t2)
t2_arr = re.split('[\n\s]+', t2_arr)[1:-1]
t2_arr = np.array(map(float, t2_arr))
t2_arr = np.reshape(t2_arr, (-1,2))
lam = t2_arr[:,0]
t2_arr = t2_arr[:,1]
t2_arr = t2_arr / np.cos(sat_zen*np.pi/180.)

# the intrinsic atmospheric reflectance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
ref_atm_arr, err3 = process.communicate(input=inp_ref_atm)
ref_atm_arr = re.split('[\n\s]+', ref_atm_arr)[1:-1]
ref_atm_arr = np.array(map(float, ref_atm_arr))
ref_atm_arr = np.reshape(ref_atm_arr, (-1,2))
lam = ref_atm_arr[:,0]
ref_atm_arr = ref_atm_arr[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# the spherical albedo
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
sph_alb_arr, err4 = process.communicate(input=inp_sph_alb)
sph_alb_arr = re.split('[\n\s]+', sph_alb_arr)[1:-1]
sph_alb_arr = np.array(map(float, sph_alb_arr))
sph_alb_arr = np.reshape(sph_alb_arr, (-1,2))
lam = sph_alb_arr[:,0]
sph_alb_arr = sph_alb_arr[:,1]

# calc of apparent TOA reflectance using separate components
app_refl = ref_atm_arr + t1_arr*t2_arr*surf_refl / (1 - sph_alb_arr*surf_refl)

# apparent TOA reflectance straight from uvspec
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
app_refl_true, err5 = process.communicate(input=inp_app_refl_true)
app_refl_true = re.split('[\n\s]+', app_refl_true)[1:-1]
app_refl_true = np.array(map(float, app_refl_true))
app_refl_true = np.reshape(app_refl_true, (-1,2))
lam = app_refl_true[:,0]
app_refl_true = app_refl_true[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# temporary plot of albedo
surf_refl = surf_refl * np.ones_like(lam)

# downsample 
lam = lam[::9]
surf_refl = surf_refl[::9]
t1_arr = t1_arr[::9]
t2_arr = t2_arr[::9]
app_refl = app_refl[::9]
ref_atm_arr = ref_atm_arr[::9]
sph_alb_arr = sph_alb_arr[::9]

plt.plot(lam, surf_refl, 'r-', label='surface reflectance')
plt.plot(lam, t1_arr*t2_arr/5, 'b--', label='transmittance/5')
plt.plot(lam, app_refl, 'g-', label='apparent reflectance')
#plt.plot(lam/1000, app_refl_true, 'r--', label='app. reflectance true')
plt.plot(lam, ref_atm_arr*2, 'b-', label='atmospheric path refl*2')
plt.plot(lam, sph_alb_arr, 'b:', label='spherical albedo')
plt.legend(loc='best')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Various (unitless)')
plt.xlim(680, wl_max)
plt.ylim(0., 0.2)
plt.show()