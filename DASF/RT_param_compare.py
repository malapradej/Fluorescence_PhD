#!/usr/bin/python
"""
Script to do comparison between 6S RT and libRadtran solutions for
transmittance, TOA reflectance, Atmospheric path reflectance, and
spherical albedo.
"""

from Py6S import *
import numpy as np
import matplotlib.pylab as plt
import subprocess
import re
import pickle
from glob import glob
from scipy.interpolate import interp1d
from netCDF4 import Dataset

plt.cla()

# the observation number in obs file and lut
obs_nr = 0

resol6s = 0.001
st_wl = 0.680
en_wl = 0.755

wls = np.arange(st_wl+resol6s, en_wl, resol6s) # remove end for interpolation

# load the GOME-2 observations
lut_files = glob('LUT_*.p')
obs_file = 'subset_remainder.p'
obs_all = pickle.load(open(obs_file))
obs = {}
for i in obs_all.keys():
    obs[i] = obs_all[i][obs_nr]

# to find rel azi need to use lut file
# index of obs in file. see that o_nr exists in file.
br = False
for index, i in enumerate(lut_files):
    for j in re.findall(r'\d+', i):
        if int(j) == obs_nr:
            br = True
            break
    if br:
        break

# load the lut and remainder file into dictionaries
lut = pickle.load(open(lut_files[index]))
# convert lists to arrays
for k in lut.keys():
    lut[k] = np.array(lut[k])

sol_zen = obs['Sol_zen']
sol_azi_6s = 180.
sol_azi_lrt = 0.
sat_zen = obs['Sat_zen']

# select observation indices in obs file
goodlut = np.where(abs(lut['sat_zen'] - sat_zen) < 0.0001)[0][0]
sat_azi = lut['rel_azi'][goodlut]
alt = lut['alt'][goodlut]
clf = obs['Clf']
st_wl_lrt = 540.
en_wl_lrt = en_wl*1000.
resollrt = 2.
wvc = 20. #40.
aot = 0.2
press = 950. #hPa
day = 25
month = 1

# load LER dataset
# need this class for script to not crash in Spyder
class DDataset():
    def __init__(self,fp):
        self.inputNC = Dataset(fp,'r')

fn = "/home/malapradej/Documents/PhD_UCL/Data/GOME-2/LER/GOME2_MetOp-A_PMD_surface_LER_product.he5"
dataset = Dataset(fn) #DDataset(fn).inputNC
lat = obs['Lat']
lon = obs['Lon']
lats = dataset.variables['Latitude'][:]
lons = dataset.variables['Longitude'][:]
idx = np.argmin(np.abs(lons - lon))
idy = np.argmin(np.abs(lats - lat))
wls_ler_pmd = dataset.variables['Wavelength'][:]
ler_min_pmd = dataset.variables['Minimum_LER'][idy,idx,:,month-1]
ler_mode_pmd = dataset.variables['Mode_LER'][idy,idx,:,month-1]

fn = "/home/malapradej/Documents/PhD_UCL/Data/GOME-2/LER/GOME2_MetOp-A_MSC_surface_LER_product.he5"
dataset = Dataset(fn) #DDataset(fn).inputNC
lat = obs['Lat']
lon = obs['Lon']
lats = dataset.variables['Latitude'][:]
lons = dataset.variables['Longitude'][:]
idx = np.argmin(np.abs(lons - lon))
idy = np.argmin(np.abs(lats - lat))
wls_ler_msc = dataset.variables['Wavelength'][:]
ler_min_msc = dataset.variables['Minimum_LER'][idy,idx,:,month-1]
ler_mode_msc = dataset.variables['Mode_LER'][idy,idx,:,month-1]

#dataset = False

lam_gome = np.array(obs['Lam_gome'])
ind_range = (lam_gome >= st_wl*1000.) * (lam_gome <= en_wl*1000.)
lam_gome = lam_gome[ind_range]
irr_gome = np.array(obs['Irr_toa'])[ind_range]
rad_gome = np.array(obs['Rad_toa'])[ind_range]
# toa refl needs to be multiplied by pi
toa_refl_gome = rad_gome / irr_gome / np.cos(sol_zen * np.pi / 180.) * np.pi

def surf_refl_fun(toa_refl, atm_refl, tr2, sph_alb):
    '''calculate boa or surface reflectance.'''
    dr = toa_refl - atm_refl
    rsurf = dr/(tr2 + sph_alb*dr)
    return rsurf

def toa_refl_fun(surf_refl, atm_refl, tr2, sph_alb):
    '''calculate toa reflectance.'''
    rtoa = atm_refl + tr2 * surf_refl / (1 - sph_alb * surf_refl)
    return rtoa

#=============================================================================
# the 6S part

s = SixS()
s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)
s.altitudes.set_target_custom_altitude(alt)
s.altitudes.set_target_pressure(press/1000.) # mb
s.altitudes.set_sensor_satellite_level()
s.geometry = Geometry.User()
s.geometry.solar_a = sol_azi_6s
s.geometry.solar_z = sol_zen
s.geometry.view_a = sat_azi
s.geometry.view_z = sat_zen
s.geometry.day = day
s.geometry.month = month
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
s.aot550 = aot
#s.ground_reflectance = GroundReflectance.HomogeneousLambertian(surf_refl)
#s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(surf_refl)
s.atmos_profile = AtmosProfile.UserWaterAndOzone(wvc*1000./10000., 0.45)
s.wavelength = Wavelength(st_wl, en_wl)

output_name='total_gaseous_transmittance'
wv, tr_gass = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(tr_gass[0])!=np.float64:
    tr_gass = np.array([tr_gass[i].total for i in np.arange(0, len(tr_gass))])

output_name='transmittance_total_scattering'
wv, tr_aero = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(tr_aero[0])!=np.float64:
    tr_aero = np.array([tr_aero[i].total for i in np.arange(0, len(tr_aero))])

# total transmittance is product of gaseous and scattering transmittances
tr2_6s = tr_gass * tr_aero

# spherical albedo (xc)
output_name='spherical_albedo'
wv, sph_alb_6s = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(sph_alb_6s[0])!=np.float64:
    sph_alb_6s = np.array([sph_alb_6s[i].total for i in np.arange(0, len(sph_alb_6s))])

# atmospheric reflectance (xb)
# 6S intrinsic atmospheric reflectance is not valid at high zenith and
# short wavelengths due to the simplification  of higher order scattering
# in 6S..... this may be the reason for differences with libRadtran....
output_name='atmospheric_intrinsic_reflectance'
wv, atm_refl_6s = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(atm_refl_6s[0])!=np.float64:
    atm_refl_6s = np.array([atm_refl_6s[i].total for i in np.arange(0, len(atm_refl_6s))])

# apparent reflectance
output_name='apparent_reflectance'
wv, toa_refl_true_6s = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(toa_refl_true_6s[0])!=np.float64:
    toa_refl_true_6s = np.array([toa_refl_true_6s[i].total for i in np.arange(0, len(toa_refl_true_6s))])

# solar irradiance
output_name='solar_spectrum'
wv, sol_irr_6s = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(sol_irr_6s[0])!=np.float64:
    sol_irr_6s = np.array([sol_irr_6s[i].total for i in np.arange(0, len(sol_irr_6s))])

toa_refl_6s = interp1d(lam_gome/1000., toa_refl_gome, kind='cubic')(wv)
surf_refl_6s = surf_refl_fun(toa_refl_6s, atm_refl_6s, tr2_6s, sph_alb_6s)
#surf_refl_6s = np.array([surf_refl]*len(wv))
#toa_refl_6s = toa_refl_fun(surf_refl_6s, atm_refl_6s, tr2_6s, sph_alb_6s)

#=============================================================================
# the libRadtran part

inp_default = 'data_files_path /usr/local/share/libRadtran/data/\n\
    output_quantity transmittance \nmol_abs_param reptran\n\
    atmosphere_file tropics \nrte_solver disort \naerosol_default \
    \naerosol_species_file continental_average \
    \naerosol_set_tau_at_wvl 550 %.3f \nmol_modify H2O %.3f MM \
    \nwavelength 400.0 1400.0 \npressure %.3f\
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran'\
    %(aot, wvc, press)

inp_t = '\nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
    \nalbedo %.1f \naltitude %.3f \nzout %s \nwavelength %.8f %.8f \
    \noutput_user lambda %s'

inp_t1 = inp_t % (sol_zen, sol_azi_lrt, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt, 'sur', st_wl_lrt, en_wl_lrt+1, 'eglo')
inp_t2 = inp_t % (sat_zen, sat_azi, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt, 'sur', st_wl_lrt, en_wl_lrt+1, 'eglo') #uu
inp_atm_refl = inp_t % (sol_zen, sol_azi_lrt, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt, 'toa', st_wl_lrt, en_wl_lrt+1, 'uu')
inp_sph_alb = inp_t % (sol_zen, sol_azi_lrt, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 1.0, alt, 'sur', st_wl_lrt, en_wl_lrt+1, 'spher_alb')
inp_sph_alb = inp_sph_alb + '\ndisort_spherical_albedo' # bug needs this at end...
#inp_toa_refl_true = inp_t % (sol_zen, sol_azi_lrt, np.cos(sat_zen*np.pi/180.)\
#    , sat_azi, surf_refl, alt, 'toa', st_wl_lrt, en_wl_lrt+1, 'uu')
inp_cloud_refl = inp_t % (sol_zen, sol_azi_lrt, np.cos(sat_zen*np.pi/180.)\
    , sat_azi, 0.0, alt, 'toa', st_wl_lrt, en_wl_lrt+1, 'uu')
inp_cloud_refl = inp_cloud_refl + '\nwc_file 1D \
    /usr/local/share/libRadtran/examples/WCSIMPLE.DAT \ncloudcover wc 1'
inp_t1 = inp_default + inp_t1
inp_t2 = inp_default + inp_t2
inp_atm_refl = inp_default + inp_atm_refl
inp_sph_alb = inp_default + inp_sph_alb
inp_cloud_refl = inp_default + inp_cloud_refl
#inp_toa_refl_true = inp_default + inp_toa_refl_true

# these are the shell pipes needed to run the command line uvspec program
# the downward transmittance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
t1, err1 = process.communicate(input=inp_t1)
t1 = re.split('[\n\s]+', t1)[1:-1]
t1 = np.array(map(float, t1))
t1 = np.reshape(t1, (-1,2))
lam = t1[:,0]
t1 = t1[:,1]
# take in consideration the cosine of solar zenith see eq(6.6) in manual
t1 = t1 / np.cos(sol_zen*np.pi/180.)

# the upward transmittance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
t2, err2 = process.communicate(input=inp_t2)
t2 = re.split('[\n\s]+', t2)[1:-1]
t2 = np.array(map(float, t2))
t2 = np.reshape(t2, (-1,2))
lam = t2[:,0]
t2 = t2[:,1]
t2 = t2 / np.cos(sat_zen*np.pi/180.)

# two way transmittance product of downward and upward transmittance
tr2_lrt = t1*t2

# the intrinsic atmospheric reflectance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
atm_refl_lrt, err3 = process.communicate(input=inp_atm_refl)
atm_refl_lrt = re.split('[\n\s]+', atm_refl_lrt)[1:-1]
atm_refl_lrt = np.array(map(float, atm_refl_lrt))
atm_refl_lrt = np.reshape(atm_refl_lrt, (-1,2))
lam = atm_refl_lrt[:,0]
atm_refl_lrt = atm_refl_lrt[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# the spherical albedo
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
sph_alb_lrt, err4 = process.communicate(input=inp_sph_alb)
sph_alb_lrt = re.split('[\n\s]+', sph_alb_lrt)[1:-1]
sph_alb_lrt = np.array(map(float, sph_alb_lrt))
sph_alb_lrt = np.reshape(sph_alb_lrt, (-1,2))
lam = sph_alb_lrt[:,0]
sph_alb_lrt = sph_alb_lrt[:,1]

# cloud reflectance
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
cld_refl_lrt, err5 = process.communicate(input=inp_cloud_refl)
cld_refl_lrt = re.split('[\n\s]+', cld_refl_lrt)[1:-1]
cld_refl_lrt = np.array(map(float, cld_refl_lrt))
cld_refl_lrt = np.reshape(cld_refl_lrt, (-1,2))
lam = cld_refl_lrt[:,0]
cld_refl_lrt = cld_refl_lrt[:,1] * np.pi / np.cos(sol_zen*np.pi/180.) * clf

# calc of apparent TOA reflectance using separate components
#toa_refl_lrt = toa_refl_fun(surf_refl, atm_refl_lrt, tr2_lrt, sph_alb_lrt)

'''
# apparent TOA reflectance straight from uvspec
process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
    subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
toa_refl_lrt_true, err5 = process.communicate(input=inp_toa_refl_true)
toa_refl_lrt_true = re.split('[\n\s]+', toa_refl_lrt_true)[1:-1]
toa_refl_lrt_true = np.array(map(float, toa_refl_lrt_true))
toa_refl_lrt_true = np.reshape(toa_refl_lrt_true, (-1,2))
lam = toa_refl_lrt_true[:,0]
toa_refl_lrt_true = toa_refl_lrt_true[:,1] * np.pi / np.cos(sol_zen*np.pi/180.)

# temporary plot of albedo
surf_refl_lrt = surf_refl * np.ones_like(lam)
'''

'''
# downsample
lam = lam[::9]
#surf_refl_lrt = surf_refl_lrt[::9]
tr2_lrt = tr2_lrt[::9]
atm_refl_lrt = atm_refl_lrt[::9]
sph_alb_lrt = sph_alb_lrt[::9]
#toa_refl_lrt_true = toa_refl_lrt_true[::9]
'''

# interpolate libRadtran data to GOME-2 wavelenghts
tr2_lrt = interp1d(lam, tr2_lrt, kind='linear')(lam_gome)
atm_refl_lrt = interp1d(lam, atm_refl_lrt, kind='linear')(lam_gome)
sph_alb_lrt = interp1d(lam, sph_alb_lrt, kind='linear')(lam_gome)
cld_toa_refl_lrt = interp1d(lam, cld_refl_lrt, kind='linear')(lam_gome)
surf_refl_lrt = surf_refl_fun(toa_refl_gome, atm_refl_lrt, tr2_lrt, sph_alb_lrt)

# correct for cloud fraction influence
cldless_toa_refl_lrt = (toa_refl_gome - cld_toa_refl_lrt) / (1 - clf)
cldless_surf_refl_lrt = surf_refl_fun(cldless_toa_refl_lrt, atm_refl_lrt,\
    tr2_lrt, sph_alb_lrt)

#=============================================================================
# plots

wv = wv*1000.

plt.plot(lam_gome, toa_refl_gome, 'm-', label='toa_refl_gome')

plt.plot(wv, tr2_6s/5., 'k--', label='tr2_6s/5')
plt.plot(wv, sph_alb_6s, 'b--', label='sph_alb_6s')
plt.plot(wv, atm_refl_6s*4., 'c--', label='atm_refl_6s*4')
#plt.plot(wv, toa_refl_6s, 'm--', label='toa_refl_6s')
plt.plot(wv, surf_refl_6s, 'r--', label='surf_refl_6s')
#plt.plot(wv, toa_refl_true_6s, 'r+', label='toa_refl_true_6s')
#plt.plot(wv, sol_irr_6s/1000., label='sol_irr_6s W/m2/nm')
plt.plot(lam_gome, tr2_lrt/5, 'k', label='tr2_lrt/5')
plt.plot(lam_gome, sph_alb_lrt, 'b', label='sph_alb_lrt')
plt.plot(lam_gome, atm_refl_lrt*4, 'c', label='atm_refl_lrt*4')
#plt.plot(lam, toa_refl_lrt, 'm', label='toa_refl_lrt')
plt.plot(lam_gome, surf_refl_lrt, 'r', label='surf_refl_lrt')
plt.plot(wls_ler_pmd, ler_min_pmd, 'g:', label='LER_PMD_MIN')
plt.plot(wls_ler_pmd, ler_mode_pmd, 'g-.', label='LER_PMD_MODE')
plt.plot(wls_ler_msc, ler_min_msc, 'm:', label='LER_MSC_MIN')
plt.plot(wls_ler_msc, ler_mode_msc, 'm-.', label='LER_MSC_MODE')
plt.plot(lam_gome, cld_toa_refl_lrt, 'm+', label='cld_toa_refl_lrt')
plt.plot(lam_gome, cldless_surf_refl_lrt, 'r+', \
    label='cldless_surf_refl_lrt')
#plt.plot(lam, toa_refl_lrt_true, 'rx', label='toa_refl_lrt_true')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Unitless ref/trans/alb')
plt.xlim(st_wl*1000., en_wl_lrt)
plt.ylim(0., 0.4)
plt.legend(loc='best', ncol=3)
plt.title('obs_nr: %d, sat_zen: %.4f, sol_zen: %.4f, clf: %.3f' % (obs_nr, \
    sat_zen, sol_zen, clf))
plt.show()
