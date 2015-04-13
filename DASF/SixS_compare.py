#!/usr/bin/python
"""
Script to do comparison with 6S RT solutions for transmittance, TOA reflectance,
Atmospheric path reflectance, and spherical albedo.
"""

from Py6S import *
import numpy as np
import matplotlib.pylab as plt

plt.cla()

resol = 0.002
st_wl = 0.680
en_wl = 0.755+resol
surf_ref = 0.1
wls = np.arange(st_wl, en_wl, resol)

s = SixS()
s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)
s.altitudes.set_target_custom_altitude(0.15)
s.altitudes.set_target_pressure(950)
s.altitudes.set_sensor_satellite_level()
s.geometry = Geometry.User()
s.geometry.solar_a = 0.
s.geometry.solar_z = 40.75
s.geometry.view_a = 342.24
s.geometry.view_z = 31.23
s.geometry.day = 25
s.geometry.month = 1
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
s.aot550 = 0.2
s.ground_reflectance = GroundReflectance.HomogeneousLambertian(surf_ref)
s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(surf_ref)
s.atmos_profile = AtmosProfile.UserWaterAndOzone(40.*1000./10000., 0.45)
s.wavelength = Wavelength(st_wl, en_wl)

'''
total transmittance downward:
T = exp(-tau/mu) + td
td = E_diff / mu / Es

where:
tau -> optical_depth
mu -> cos(geometry.solar_z)
E_diff -> diffuse_solar_irradiance
Es -> solar_spectrum

total transmittance upward:
T = exp(-tau/mu) + td
'''

# total_gaseous_transmittance
# transmittance_aerosol_scattering
# transmittance_total_scattering

output_name='total_gaseous_transmittance'
wv, tr_gass = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(tr_gass[0])!=np.float64:
    tr_gass = np.array([tr_gass[i].total for i in np.arange(0, len(tr_gass))])
#SixSHelpers.Wavelengths.plot_wavelengths(wv, tr2, output_name)

output_name='transmittance_total_scattering'
wv, tr_aero = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(tr_aero[0])!=np.float64:
    tr_aero = np.array([tr_aero[i].total for i in np.arange(0, len(tr_aero))])
tr2 = tr_gass * tr_aero
#SixSHelpers.Wavelengths.plot_wavelengths(wv, tr2, output_name)

# spherical albedo (xc)
output_name='spherical_albedo'
wv, sph_alb = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(sph_alb[0])!=np.float64:
    sph_alb = np.array([sph_alb[i].total for i in np.arange(0, len(sph_alb))])
#SixSHelpers.Wavelengths.plot_wavelengths(wv, sph_alb, output_name)

# atmospheric reflectance (xb)
# the intrinsic atmospheric reflectance is not valid at high zenith and
# short wavelengths due to the simplification  of higher order scattering
# in 6S..... this is reason for differences with libRadtran....
output_name='atmospheric_intrinsic_reflectance'
wv, atm_int_ref = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(atm_int_ref[0])!=np.float64:
    atm_int_ref = np.array([atm_int_ref[i].total for i in np.arange(0, len(atm_int_ref))])
#SixSHelpers.Wavelengths.plot_wavelengths(wv, atm_int_ref, output_name)

# apparent reflectance
output_name='apparent_reflectance'
wv, app_ref_true = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(app_ref_true[0])!=np.float64:
    app_ref_true = np.array([app_ref_true[i].total for i in np.arange(0, len(app_ref_true))])

# solar irradiance
output_name='solar_spectrum'
wv, sol_irr = SixSHelpers.Wavelengths.run_wavelengths(s, wls, output_name=output_name)
if type(sol_irr[0])!=np.float64:
    sol_ir = np.array([sol_irr[i].total for i in np.arange(0, len(sol_irr))])
#app_ref = atm_int_ref + tr2*surf_ref / (1. - sph_alb*surf_ref)

def surf_refl(app_ref_true, atm_int_ref, tr2, sph_alb):
    '''calculate boa or surface reflectance.'''
    dr = app_ref_true - atm_int_ref
    rsurf = dr/(tr2 + sph_alb*dr)
    return rsurf

def toa_refl(surf_ref, atm_int_ref, tr2, sph_alb):
    '''calculate toa reflectance.'''
    rtoa = atm_int_ref + tr2 * surf_ref / (1 - sph_alb * surf_ref)
    return rtoa

#surf_ref = surf_refl(app_ref_true, atm_int_ref, tr2, sph_alb)
surf_ref = np.array([surf_ref]*len(wv))
toa_ref = toa_refl(surf_ref, atm_int_ref, tr2, sph_alb)

plt.plot(wv, tr2/5., 'b--', label='dbl_trans/5')
plt.plot(wv, sph_alb, 'b:', label='sph_alb')
plt.plot(wv, atm_int_ref*2., 'b-', label='atm_path*2')
plt.plot(wv, toa_ref, 'g-', label='toa_ref')
plt.plot(wv, surf_ref, 'r-', label='surf_ref')
#plt.plot(wv, sol_irr/1000., label='sol_irr W/m2/nm')
plt.legend(loc='best')
plt.ylabel('Unitless ref/trans/alb')
plt.xlabel('Wavelength (microns)')
plt.ylim(0., 0.2)
plt.xlim(st_wl, en_wl-resol)
plt.show()
#SixSHelpers.Wavelengths.plot_wavelengths(wv, app_ref, 'apparent_reflectance')
#SixSHelpers.Wavelengths.plot_wavelengths(wv, app_ref_true, 'apparent_reflectance_true')
