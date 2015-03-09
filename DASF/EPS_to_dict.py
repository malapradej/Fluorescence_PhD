#!/usr/bin/python
"""
A script that creates a pickles dict file from GOME-2 l1B data. Requires one
argument which is the path and filename of the EPS format GOME-2 l1b file.
Usage: $ ./DASF_GOME-2.py </path/to/file>
The script creates a pickled dictionary file per orbit filtered according to
the QA limits set out below. 

The script requires CODA library for python to be installed.

To run the script certain directories need to be set up in advance. They include
directories created manualy or during installation of CODA library. The following
directories need to be created and amended in the code below:
variable            usual path can be amended
========            =========================
path_coda           /usr/local/share/coda/definitions/
path_dasf           path for output *.p pickled dictionary DASF files.
"""
import sys
import os
# sets the CODA definition file path. This is used by CODA to understand the
# format of the EPS file.
path_coda = '/usr/local/share/coda/definitions/'
if not(os.path.isdir(path_coda)):
    text = 'The path %s does not exist.' % path_coda
    raise IOError(text)
os.putenv('CODA_DEFINITION', path_coda)
import coda
#import beatl2
import numpy as np
#import matplotlib.pylab as plt
import time
#import subprocess
#from scipy.interpolate import interp1d
#from scipy import stats
#import matplotlib as mpl
import pickle

def mmol_Watt(mmols, wl):
    '''A function that converts a measurement in micro mol photons to Watts. 
    This is usefull in flux calculations. 
    Input: mmols - mmols at wavelength, wl - wavelength in nm.
    Output: Watts    
    The conversion is based on the examples in 
    http://5e.plantphys.net/article.php?ch=t&id=131
    '''
    E = 1.988e-16 / wl
    micro = 10.**-6
    avo = 6.02e23
    # convert micromoles (umol) per m2 per sec to W per m2
    Watt = np.multiply(mmols*micro*avo, E)
    return Watt
    
def phot_Watt(photons, wl):
    '''Function to convert radiant flux from photons / (cm2 s nm) to 
    Watt / m2.
    Input: photons - photons at wavelenght, wl - wavelength in nm.
    Output: Watt - Watt / m2
    '''
    # convert to m^2 from of cm^2
    photons *= 10000.
    # convert the solar irradiance from photons to mols
    avo = 6.02e23
    mols = photons / avo 
    # convert to mmols
    mmols = mols * 1.0e6
    Watt = mmol_Watt(mmols, wl)
    return Watt

def interpolate(rc, ic):
    '''Function that interpolates a solar irradiance spectrum to the same
    wavelengths as the radiance spectrum. This is due to slight differences
    between the radiance and irradiance wavelength intervals.
    This function does not get passed the radiance and irradiance arrays but
    only the indexes of the radiance wavelenght at which to interpolate, and 
    the index to the upper irradiance wavelength used in the interpolation.
    Input: rc - radiance index, ic - irradiance index.
    Output: inter_irr - interpolated irradiance value.
    '''
    lam_int = sol_lam[ic] - sol_lam[ic-1]
    step = rad_lam[rc] - sol_lam[ic-1]
    frac = step / lam_int
    irr_int = sol_irr[ic] - sol_irr[ic-1]
    inter_irr = sol_irr[ic-1] + frac * irr_int
    return inter_irr
    
# range of forward scans out of [0:31]
forw_scans = np.arange(24)

# cloud fraction threshold
cloud_thresh = 0.1

# solar zenith angle threshold
sol_zen_lim = 80.0

# satelite zenith angle threshold
sat_zen_lim = 70.0

# path and EPS file to read
if len(sys.argv)!=2:
    text = '1 /path/file argument required and %d given.' % len(sys.argv)-1
    raise IOError(text)
e_fn = sys.argv[1]
if not(os.path.exists(e_fn)):
    text = 'The file %s does not exist.' % e_fn
    raise IOError(text)
#e_fn = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/l1b/2007/01/GOME_xxx_1B_M02_20070125124939Z_20070125143057Z_R_O_20120209112017Z'
# opening the file
ef = coda.open(e_fn)

# creates a base name without extension used in rest of code
base_fn = os.path.splitext(e_fn)[0] # removes any extension from path
base_fn = os.path.basename(base_fn) # selects only filename

# path to pickles dict data to be saved
path_dasf = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/DASF/'
if not(os.path.exists(path_dasf)):
    text = 'The path %s does not exist.' % path_dasf
    raise IOError(text)

# reading solar wavelength and irradiance in band 4
sol_lam = coda.fetch(ef, 'VIADR_SMR', 0, 'LAMBDA_SMR')[3]
n_sol_lam = len(sol_lam)
sol_irr = coda.fetch(ef, 'VIADR_SMR', 0, 'SMR')[3]
sol_irr = phot_Watt(sol_irr, sol_lam)

# reading total number of scans
n_mdr = coda.fetch(ef, 'MPHR', 'TOTAL_MDR')

# list for point lat, lon, radiance, irradiance, spectrum etc.
Lat = []
Lon = []
Spec = []
Rad = []
Irr = []
Ref_toa = []
Time = []
Alt = []
Cloud_frac = []
Sol_zen = []
Sol_azi = []
Sat_zen = []
Sat_azi = []

for cur_scan in np.arange(n_mdr):
    sys.stdout.write('\rcurrent scan line is %d out of %d scans.' % \
        (cur_scan+1, n_mdr))
    sys.stdout.flush()
    # is this an Earthshine scan, not calibration etc, then skip this scan?
    es = coda.get_field_available(ef, 'MDR', cur_scan, 'Earthshine')
    if not(es):
        continue
    
    # is this a nadir scan (0), not calibration etc?
    nad = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'OBSERVATION_MODE')
    nad = True if nad == 0 else False
    
    # QA of scan.
    # is scan degraded due to instrument degradation?
    bad_ins = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'DEGRADED_INSTR_MDR')
    bad_ins = True if bad_ins == 1 else False
    # is scan degraded due to processing degradation?
    bad_pro = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'DEGRADED_PROC_MDR')
    bad_pro = True if bad_pro == 1 else False
    # if some processing degraded data is this critical for our purpose?
    if bad_pro:
        # is it due to missing data which is not ok, or old calibration data
        # which is still ok?
        f_miss = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'PCD_BASIC',\
            'F_MISS')
        bad_pro = True if f_miss != 0 else False
        
    
    # is this calibrated radiance (0) or reflectance (1)?
    rad_ref = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'OUTPUT_SELECTION')
    rad_ref = True if rad_ref == 0 else False
    
    # if any above conditions are not met skip this scan
    if not(es) or not(nad) or bad_ins or bad_pro or not(rad_ref):
        continue
        
    for cur_pnt in forw_scans:  
        # did cloud fraction retrieval work, if not skip point?
        cf = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'CLOUD', 'FAIL_FLAG',\
            cur_pnt)
        cf = True if cf == 0 else False
        if not(cf):
            continue
        
        # is parameter for cloud fraction (0) not snow/ice (1), if not skip point?
        cf = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'CLOUD', 'FIT_MODE',\
            cur_pnt)
        cf = True if cf == 0 else False
        if not(cf):
            continue
        
        # reading cloud fraction, and if > threshold skip point
        cf = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'CLOUD', 'FIT_2',\
            cur_pnt)
        if cf > cloud_thresh:
            continue
        
        # reading solar zenith angle if greater than threshold skip point
        sol_zen = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'SOLAR_ZENITH', [1, cur_pnt])
        if sol_zen > sol_zen_lim:
            continue
        
        # satellite zenith angle
        sat_zen = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'SAT_ZENITH', [1, cur_pnt])
        if sat_zen > sat_zen_lim:
            continue
        
        # number of spectral records in band 4, and create point array for radiance
        n_spec = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'REC_LENGTH', 5)
        rad_arr = np.zeros(n_spec)
        
        # placeholder for dud spectrum. Check against this after loop
        dud = False
        
        # reading coordinates
        lat, lon  = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'CENTRE', cur_pnt)
        
        
        # TEMPORARY halt to see Amazon....
        '''
        if lat < 5.0 and lon < -53.0:
            pass
        else:
            continue
        '''
        
        for cur_spc in np.arange(n_spec):
            # reading radiance, error and stokes fraction
            rad, rad_e, st_f  = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'BAND_4',\
                [cur_pnt,cur_spc])
            
            # is the current point radiance nan or negative, ie ocean etc, then
            # skip the whole point?
            dud = True if rad < 0. or np.isnan(rad) else False
            if dud:
                break
            
            rad_arr[cur_spc] = rad
            
        # skip point if dud in radiance array
        if dud:
            continue
        
        # satellite azimuth angle
        sat_azi = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'SAT_AZIMUTH', [1, cur_pnt])
        
        # reading solar azimuth angle
        sol_azi = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'SOLAR_AZIMUTH', [1, cur_pnt])
        
        # reading in time of observation
        time_s = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_BASIC',\
            'UTC_TIME', cur_pnt)
        # convert time from 2000-01-01 to 1970-01-01 epoch
        date_time = '2000-01-01 00:00:00'
        pattern = '%Y-%m-%d %H:%M:%S'
        offset = time.mktime(time.strptime(date_time, pattern))
        time_s += offset
        
        # reading in altitude of point
        alt = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'GEO_EARTH',\
            'SURFACE_ELEVATION', cur_pnt)
        
        # reading wavelength for radiance data
        rad_lam = coda.fetch(ef, 'MDR', cur_scan, 'Earthshine', 'WAVELENGTH_4')

        # convert radiance and radiance error from photons to Watts
        rad_arr = phot_Watt(rad_arr, rad_lam)
        
        # interpolate the solar irradiance data to the same wavelengths as the
        # radiance data
        n_rad_lam = len(rad_lam)
        n = 0
        ic = 0
        irr = []
        while rad_lam[n] < sol_lam[0]:
            n += 1
        for rc in np.arange(n, n_rad_lam):
            while sol_lam[ic] < rad_lam[rc] and ic != n_sol_lam-1:
                ic += 1
            if ic != n_sol_lam-1:
                irr.append(interpolate(rc, ic))
            else:
                if rc != n_rad_lam-1:
                    if rad_lam[rc+1] < sol_lam[ic]:
                        irr.append(interpolate(rc, ic))
                        continue
                    else:
                        irr.append(interpolate(rc, ic))
                        break
                else:
                    irr.append(interpolate(rc, ic))
                    break
        # copy back only the values corresponding to interpolated wavelengths
        # into radiance, coordinate and wavelength arrays
        rad_arr = rad_arr[n: rc+1]
        rad_lam = rad_lam[n: rc+1]
        
        # calculate reflectance at TOA. 
        irr_arr = np.array(irr)
        ref_toa_arr = rad_arr / (irr_arr * np.cos(sol_zen*np.pi/180.))
        
        # append point data to orbit lists
        #DASF.append(dasf)
        #QA_dasf.append((slope, intercept, r_value, p_value))
        Lat.append(lat)
        Lon.append(lon)
        Rad.append(rad_arr)         # radiance at TOA
        Spec.append(rad_lam)        # wavelenghts of the radiance
        Irr.append(irr_arr)         # irradiance at TOA
        Ref_toa.append(ref_toa_arr) # reflectance at TOA
        Time.append(time_s)         # time in sec since 1970-01-01
        Alt.append(alt)             # altitude in meters
        Cloud_frac.append(cf)       # cloud fraction
        Sol_zen.append(sol_zen)
        Sol_azi.append(sol_azi)
        Sat_zen.append(sat_zen)
        Sat_azi.append(sat_azi)

# new line after stdout
sys.stdout.write("\n")        
# closes file to free memory
coda.close(ef)

# saves the data into a dictionary on disk
disk_dict = dict(Lat=Lat, Lon=Lon, Rad=Rad,\
    Spec=Spec, Irr=Irr, Ref_toa=Ref_toa, Time=Time,\
    Alt=Alt, Cloud_frac=Cloud_frac, Sol_zen=Sol_zen, Sol_azi=Sol_azi,\
    Sat_zen=Sat_zen, Sat_azi=Sat_azi)
d_fn = os.path.join(path_dasf, base_fn + '.p')
pf = open(d_fn, 'wb')
pickle.dump(disk_dict, pf)
pf.close()

print 'total of #d points saved to disk.' %len(Lat)

# unloads the CODA module
coda.done()