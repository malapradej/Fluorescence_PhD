#!/usr/bin/python
"""
A script to calculate libRadtran transmittances, spherical albedo, intrinsic atmopsheric
reflectance and apparent reflectance to compare with 6S.
"""

import subprocess
import matplotlib.pylab as plt
import numpy as np
import re

sol_zen = 60.
sol_azi = 90.
sat_zen = 60.
sat_azi = 270.
alt = 0.
wl_min = 400.
wl_max = 800.

surf_alb = 0.5

inp_default = 'data_files_path /usr/local/share/libRadtran/data/\n\
    output_quantity transmittance \nmol_abs_param reptran\n\
    atmosphere_file US-standard \nrte_solver disort \naerosol_default \
    \naerosol_species_file continental_average \
    \naerosol_set_tau_at_wvl 550 0.12 \nmol_modify H2O 9.7 MM \
    \nmol_modify O3 450 DU \nwavelength 400.0 1400.0 \
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran'

t_inp = '\nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
    \nalbedo %.1f \naltitude %.3f \nzout %s \nwavelength %.8f %.8f \
    \noutput_user lambda %s'
    
t1_inp = t_inp % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'sur', wl_min, wl_max, 'eglo')
t2_inp = t_inp % (sat_zen, sat_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'sur', wl_min, wl_max, 'eglo') #uu
ref_atm_inp = t_inp % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt/1000., 'toa', wl_min, wl_max, 'uu')
sph_alb_inp = t_inp % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 1.0, alt/1000., 'sur', wl_min, wl_max, 'spher_alb')
sph_alb_inp = sph_alb_inp + '\ndisort_spherical_albedo' # bug needs this at end...
app_refl_true_inp = t_inp % (sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, surf_alb, alt/1000., 'toa', wl_min, wl_max, 'uu')
t1_inp = inp_default + t1_inp
t2_inp = inp_default + t2_inp
ref_atm_inp = inp_default + ref_atm_inp
sph_alb_inp = inp_default + sph_alb_inp
app_refl_true_inp = inp_default + app_refl_true_inp

# these are the shell pipes needed to run the command line uvspec program
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
t2_arr = t2_arr / np.cos(sat_zen*np.pi/180.)

# the intrinsic atmospheric reflectance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
ref_atm_arr, err3 = process.communicate(input=ref_atm_inp)
ref_atm_arr = re.split('[\n\s]+', ref_atm_arr)[1:-1]
ref_atm_arr = np.array(map(float, ref_atm_arr))
ref_atm_arr = np.reshape(ref_atm_arr, (-1,2))
lam = ref_atm_arr[:,0]
ref_atm_arr = ref_atm_arr[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# the spherical albedo
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
sph_alb_arr, err4 = process.communicate(input=sph_alb_inp)
sph_alb_arr = re.split('[\n\s]+', sph_alb_arr)[1:-1]
sph_alb_arr = np.array(map(float, sph_alb_arr))
sph_alb_arr = np.reshape(sph_alb_arr, (-1,2))
lam = sph_alb_arr[:,0]
sph_alb_arr = sph_alb_arr[:,1]

# calc of apparent TOA reflectance using separate components
app_refl = ref_atm_arr + t1_arr*t2_arr*surf_alb / (1 - sph_alb_arr*surf_alb)

# apparent TOA reflectance straight from uvspec
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
app_refl_true, err5 = process.communicate(input=app_refl_true_inp)
app_refl_true = re.split('[\n\s]+', app_refl_true)[1:-1]
app_refl_true = np.array(map(float, app_refl_true))
app_refl_true = np.reshape(app_refl_true, (-1,2))
lam = app_refl_true[:,0]
app_refl_true = app_refl_true[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# temporary plot of albedo
surf_alb = surf_alb * np.ones_like(lam)
plt.plot(lam, surf_alb, 'r-', label='surface albedo')
plt.plot(lam, t1_arr*t2_arr, 'b--', label='transmittance')
plt.plot(lam, app_refl, 'g-', label='apparent reflectance')
plt.plot(lam/1000, app_refl_true, 'r--', label='app. reflectance true')
plt.plot(lam, ref_atm_arr, 'b-', label='atmospheric albedo')
plt.plot(lam, sph_alb_arr, 'b:', label='spherical albedo')
plt.legend(loc='upper left')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Various (unitless)')
plt.show()