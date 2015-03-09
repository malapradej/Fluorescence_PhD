#!/usr/bin/python
"""
Script to do comparison with 6S RT solutions for transmittance, TOA reflectance,
Atmospheric path reflectance, and spherical albedo.
"""

from Py6S import *
import numpy as np
#import matplotlib.pylab as plt

st_wl = 0.530
en_wl = 0.780
resol = 0.005
surf_alb = 0.5
wls = np.arange(st_wl, en_wl, resol)

s = SixS()
s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.USStandard1962)
s.wavelength = Wavelength(st_wl, en_wl)
s.altitudes.set_target_sea_level()
s.altitudes.set_sensor_satellite_level()
s.geometry = Geometry.User()
s.geometry.solar_a = 90
s.geometry.solar_z = 60
s.geometry.view_a = 270
s.geometry.view_z = 60
s.geometry.day = 25
s.geometry.month = 1
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
s.aot550 = 0.1
s.ground_reflectance = GroundReflectance.HomogeneousLambertian(surf_alb)
s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(surf_alb)

s.aot550 = 0.12
s.atmos_profile = AtmosProfile.UserWaterAndOzone(0.97, 0.45)

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
wv, tr_gass = SixSHelpers.Wavelengths.run_vnir(s, output_name=output_name)
if type(tr_gass[0])!=np.float64:
    tr_gass = [tr_gass[i].total for i in np.arange(0, len(tr_gass))]
tr_gass = np.array(tr_gass)
#SixSHelpers.Wavelengths.plot_wavelengths(wv, tr2, output_name)

output_name='transmittance_total_scattering'
wv, tr_aero = SixSHelpers.Wavelengths.run_vnir(s, output_name=output_name)
if type(tr_aero[0])!=np.float64:
    tr_aero = [tr_aero[i].total for i in np.arange(0, len(tr_aero))]
tr_aero = np.array(tr_aero)
tr2 = tr_gass * tr_aero
#SixSHelpers.Wavelengths.plot_wavelengths(wv, tr2, output_name)

# spherical albedo (xc)
output_name='spherical_albedo'
wv, sph_alb = SixSHelpers.Wavelengths.run_vnir(s, output_name=output_name)
if type(sph_alb[0])!=np.float64:
    sph_alb = [sph_alb[i].total for i in np.arange(0, len(sph_alb))]
sph_alb = np.array(sph_alb)
#SixSHelpers.Wavelengths.plot_wavelengths(wv, sph_alb, output_name)

# atmospheric reflectance (xb)
# the intrinsic atmospheric reflectance is not valid at high zenith and
# short wavelengths due to the simplification  of higher order scattering
# in 6S..... this is reason for differences with libRadtran....
output_name='atmospheric_intrinsic_reflectance'
wv, atm_int_ref = SixSHelpers.Wavelengths.run_vnir(s, output_name=output_name)
if type(atm_int_ref[0])!=np.float64:
    atm_int_ref = [atm_int_ref[i].total for i in np.arange(0, len(atm_int_ref))]
#SixSHelpers.Wavelengths.plot_wavelengths(wv, atm_int_ref, output_name)

# apparent reflectance
output_name='apparent_reflectance'
wv, app_ref_true = SixSHelpers.Wavelengths.run_vnir(s, output_name=output_name)
if type(app_ref_true[0])!=np.float64:
    app_ref_true = [app_ref_true[i].total for i in np.arange(0, len(app_ref_true))]

app_ref = atm_int_ref + tr2*surf_alb / (1. - sph_alb*surf_alb)

SixSHelpers.Wavelengths.plot_wavelengths(wv, app_ref, 'apparent_reflectance')
SixSHelpers.Wavelengths.plot_wavelengths(wv, app_ref_true, 'apparent_reflectance_true')
