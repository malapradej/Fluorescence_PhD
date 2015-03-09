#!/usr/bin/python
"""
This version uses wavelenghts between 690 and 710nm only, and it does not
remove atmospheric effects....
A script that calculates DASF from GOME-2 l1B data. The script requires one
argument which is the path and filename of the EPS format GOME-2 l1b file.
Usage: $ ./DASF_GOME-2.py </path/to/file>
The script creates a pickled dictionary file per orbit with all DASF values, 
observation data and QA data. The DASF calculations require a typical leaf 
spectrum created using PROSPECT in the same directiory as the script named 
leaf_spectrum.txt. It also requires the default.inp file in the same directory.

The script requires CODA library for python to be installed as well as 
libRadtran's uvspec. 

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
import matplotlib.pylab as plt
import time
import subprocess
from scipy.interpolate import interp1d
from scipy import stats
import matplotlib as mpl
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
sol_zen_lim = 70.0

# satelite zenith angle threshold
sat_zen_lim = 70.0

# lower and upper limits of the red-edge
edge_lower = 690.
# TRY to remove the water vapour absorption lines starting at about 710nm
edge_upper = 710. #730.

# default libRadtran input file text to be added to
inp_file = open('default.inp', 'r')
inp_default = inp_file.read()
inp_file.close()

# path and EPS file to read
if len(sys.argv)!=2:
    text = '1 /path/file argument required and %d given.' % len(sys.argv)-1
    raise IOError(text)
e_fn = sys.argv[1]
if not(os.path.exists(e_fn)):
    text = 'The file %s does not exist.' % e_fn
    raise IOError(text)
e_fn = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/l1b/2007/01/GOME_xxx_1B_M02_20070125124939Z_20070125143057Z_R_O_20120209112017Z'
# opening the file
ef = coda.open(e_fn)

# creates a base name without extension used in rest of code
base_fn = os.path.splitext(e_fn)[0] # removes any extension from path
base_fn = os.path.basename(base_fn) # selects only filename

# path to DASF data to be saved
path_dasf = '/home/malapradej/Documents/PhD_UCL/Data/GOME-2/DASF'
if not(os.path.exists(path_dasf)):
    text = 'The path %s does not exist.' % path_dasf
    raise IOError(text)

# leaf albedo file generated from PROSPECT 
# typical leaf level parameters generated based on metadata-analysis
# see spreadsheet in PROSPECT folder
l_lower = 580. # limits for plotting purposes
l_upper = 800.
l_fn = 'leaf_spectrum.txt'
leaf_alb_arr = np.genfromtxt(l_fn)
leaf_alb_arr = leaf_alb_arr[np.intersect1d(np.array(np.where(leaf_alb_arr[:,0]\
    >=l_lower)),np.array(np.where(leaf_alb_arr[:,0]<=l_upper)))]
leaf_lam = leaf_alb_arr[:,0]
leaf_alb_arr = leaf_alb_arr[:,1] + leaf_alb_arr[:,2]

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
Ref_sur = []
Ref_toa = []
Time = []
Alt = []
Cloud_frac = []
Sol_zen = []
Sol_azi = []
Sat_zen = []
Sat_azi = []
DASF = []
QA_dasf = []

for cur_scan in np.arange(n_mdr):
    print 'current scan line is %d out of %d scans.' % (cur_scan+1, n_mdr)
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
        # this can be replaced by ground level reflectance taking in
        # consideration atmospheric RT....
        irr_arr = np.array(irr)
        # needs to be understood as reflectance not albedo. rename later....
        ref_toa_arr = rad_arr / (irr_arr * np.cos(sol_zen*np.pi/180.))
        
        # this section is to derive surface level albedo as opposed to
        # just using the TOA values. A work in progress...
        
        
        # save the wavelengths and irradiance array to a temp disk file
        # as irradiance file for uvspec
        lam_irr_fn = 'lam_irr_%s.dat' % base_fn
        np.savetxt(lam_irr_fn, np.append([rad_lam],[irr_arr], axis=0).T, \
            fmt='%.8f')
        '''
        # make this version not consider the atmosphere ie. use TOA albedo
        # UNCOMMENT below otherwise....
        
        # calculate the transmittances down and up and atmospheric
        # albedo and spherical albedo using uvspec
        # first set up the .inp file for t1 and t2, ref_atm_arr, sph_alb_arr
        pattern = '%Y %m %d %H %M %S'
        time_txt = time.strftime(pattern, time.gmtime(time_s))
        wl_min = rad_lam[0] #-0.01
        wl_max = rad_lam[-1] #+0.01
        t_inp = '\ntime %s \nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
            \nalbedo %.1f \naltitude %.3f \nzout %s \nsource solar %s \
            \nwavelength %.8f %.8f \noutput_user %s \nlatitude %s %.6f \
            \nlongitude %s %.6f'
        ns = 'N' if lat >= 0.0 else 'S'
        ew = 'E' if lon >= 0.0 else 'W'
        t1_inp = t_inp % (time_txt, sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
            , sat_azi, 0.0, alt/1000., 'sur', lam_irr_fn, wl_min, wl_max, \
            'eglo', ns, abs(lat), ew, abs(lon))
        t2_inp = t_inp % (time_txt, sol_zen, sol_azi, np.cos(sat_zen*np.pi/180.)\
            , sat_azi, 1.0, alt/1000., 'toa', lam_irr_fn, wl_min, wl_max, \
            'uu', ns, abs(lat), ew, abs(lon))
        ref_atm_inp = t_inp % (time_txt, sol_zen, sol_azi, \
            np.cos(sat_zen*np.pi/180.), sat_azi, 0.0, alt/1000., 'toa', \
            lam_irr_fn, wl_min, wl_max, 'uu', ns, abs(lat), ew, abs(lon))
        sph_alb_inp = t_inp % (time_txt, sol_zen, sol_azi, \
            np.cos(sat_zen*np.pi/180.), sat_azi, 1.0, alt/1000., 'sur', \
            lam_irr_fn, wl_min, wl_max, 'spher_alb', ns, abs(lat), ew, abs(lon))
        sph_alb_inp = sph_alb_inp + '\ndisort_spherical_albedo' # bug needs this at end...
        t1_inp = inp_default + t1_inp
        t2_inp = inp_default + t2_inp
        ref_atm_inp = inp_default + ref_atm_inp
        sph_alb_inp = inp_default + sph_alb_inp
        
        # these are the shell pipes needed to run the command line uvspec program
        process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
            subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        t1_arr, err1 = process.communicate(input=t1_inp)
        t1_arr = t1_arr.split('\n')[:-1]
        t1_arr = np.array(map(float, t1_arr))
        # take in consideration the cosine of solar zenith see eq(6.6) in manual
        t1_arr = t1_arr / np.cos(sol_zen*np.pi/180.)
        process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
            subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        t2_arr, err2 = process.communicate(input=t2_inp)
        t2_arr = t2_arr.split('\n')[:-1]
        t2_arr = np.array(map(float, t2_arr))
        # t2 is the tau in Chandrasekhar (1960) equation (191) p.271
        # we convert it to transmittance using the exponential similar to 
        # Beers law.
        t2_arr = np.exp(-t2_arr/np.cos(sat_zen*np.pi/180.)) 
        process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
            subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        ref_atm_arr, err3 = process.communicate(input=ref_atm_inp)
        ref_atm_arr = ref_atm_arr.split('\n')[:-1]
        ref_atm_arr = np.array(map(float, ref_atm_arr))
        ref_atm_arr = ref_atm_arr * np.pi # need to multiply by pi
        process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
            subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        sph_alb_arr, err4 = process.communicate(input=sph_alb_inp)
        sph_alb_arr = sph_alb_arr.split('\n')[:-1]
        sph_alb_arr = np.array(map(float, sph_alb_arr))
        # convert to hemispherical and take in consideration the cosine of sol_zen
        # see notes on 9/01/15 eq(3)
        ref_atm_arr = ref_atm_arr * np.pi / (np.cos(sol_zen*np.pi/180.) \
            * np.cos(sat_zen*np.pi/180.))
        # ground level albedo see Tilstra & Stammes (2014) or note on 9/01/15
        ref_sur_arr = (ref_toa_arr - ref_atm_arr) / (t1_arr * t2_arr - sph_alb_arr*\
            ref_atm_arr + sph_alb_arr*ref_toa_arr)


        # temporary plot of albedo
        plt.plot(rad_lam, ref_sur_arr, 'r-', label='surface reflectance')
        plt.plot(rad_lam, t1_arr*t2_arr, 'b--', label='transmittance')
        plt.plot(rad_lam, ref_toa_arr, 'g-', label='TOA reflectance')
        plt.plot(rad_lam, ref_atm_arr, 'b-', label='atmospheric reflectance')
        plt.plot(rad_lam, sph_alb_arr, 'b:', label='spherical albedo')
        plt.plot(leaf_lam, leaf_alb_arr, 'y-', label='leaf albedo')
        plt.legend(loc='center left')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Various (unitless)')
        plt.show()
        
        
        # COMMENT/REMOVE below for with atmosphere....
        #ref_sur_arr = ref_toa_arr
        
        # calculation of DASF is done here
        # start with selecting only wavelenghts over the red-edge 
        
        index_alb = np.intersect1d(np.array(np.where(rad_lam>edge_lower)),\
            np.array(np.where(rad_lam<edge_upper)))
        index_leaf = np.intersect1d(np.array(np.where(leaf_lam>=edge_lower)),\
            np.array(np.where(leaf_lam<=edge_upper)))
        edge_alb_lam = rad_lam[index_alb]
        edge_alb_arr = alb_sur_arr[index_alb]
        edge_leaf_lam = leaf_lam[index_leaf]
        edge_leaf_arr = leaf_alb_arr[index_leaf]
        
        #interpolate leaf albedo wavelenths to surface albedo wavelengths
        fun_interp = interp1d(edge_leaf_lam, edge_leaf_arr, kind='linear')
        edge_leaf_arr = fun_interp(edge_alb_lam)
        
        # setup fitting data as described in Disney & Lewis (2005)
        y = edge_leaf_arr / edge_alb_arr 
        x = edge_leaf_arr
        
        # fit and retrieve the slope and intercept using linear method
        slope, intercept, r_value, p_value, std_err = stats.linregress(-x, y)
        
        # retrieve DASF based on Huang et al (2007) see notes on 12/12/14
        dasf = 1. / (intercept - slope)
        
        
        # plot some informative graphs
        # leaf vs surface albedo
        # COMMENT below....
        
        plt.plot(edge_leaf_arr, edge_alb_arr, 'r.')
        plt.xlabel('Leaf albedo')
        plt.ylabel('Surface albedo')
        plt.show()
        
        # leaf albedo vs leaf over surface albedo
        plt.plot(x, y, 'r.')
        plt.xlabel('Leaf albedo')
        plt.ylabel('Leaf albedo / Surface albedo')
        plt.plot(x, -x*slope + intercept, 'k--')
        text = 'slope=%.4f, intercept=%.4f\nR=%.4f, p=%.4e\nDASF=%.4f' % \
            (-slope, intercept, r_value, p_value, dasf)
        plt.text(0.5, 1.5, text)
        plt.show()
        '''
        
        # append point data to orbit lists
        #DASF.append(dasf)
        #QA_dasf.append((slope, intercept, r_value, p_value))
        Lat.append(lat)
        Lon.append(lon)
        Rad.append(rad_arr)         # radiance at TOA
        Spec.append(rad_lam)        # wavelenghts of the radiance
        Irr.append(irr_arr)         # irradiance at TOA
        Ref_toa.append(ref_toa_arr) # reflectance at TOA
        #Ref_sur.append(alb_sur_arr) # reflectance at surface
        Time.append(time_s)         # time in sec since 1970-01-01
        Alt.append(alt)             # altitude in meters
        Cloud_frac.append(cf)       # cloud fraction
        Sol_zen.append(sol_zen)
        Sol_azi.append(sol_azi)
        Sat_zen.append(sat_zen)
        Sat_azi.append(sat_azi)
        
# closes file to free memory
coda.close(ef)

# plots the DASF values
'''
plt.scatter(Lon, Lat, c=DASF, s=50, cmap=mpl.cm.gray)
plt.show()
'''

# saves the data into a dictionary on disk
disk_dict = dict(Lat=Lat, Lon=Lon, Rad=Rad,\
    Spec=Spec, Irr=Irr, Ref_toa=Ref_toa, Time=Time,\
    Alt=Alt, Cloud_frac=Cloud_frac, Sol_zen=Sol_zen, Sol_azi=Sol_azi,\
    Sat_zen=Sat_zen, Sat_azi=Sat_azi)
d_fn = os.path.join(path_dasf, base_fn + '.p')
pf = open(d_fn, 'wb')
pickle.dump(disk_dict, pf)
pf.close()

# removes the temporary irradiance file on disk
os.remove(lam_irr_fn)

# unloads the CODA module
coda.done()