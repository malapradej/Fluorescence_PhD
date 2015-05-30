#!/usr/bin/python
"""
A script to do PCA of cloud reflectance data. PCA will be used in the retrieval
of other parameters.
A PCA will cover each view and solar zenith as per angles in LUT file.
The samples (features) used to cover the parameter space will be selected based
on a literature review of global frequency of types as well as properties. The
types frequency and properties are found in the file PCA_cloud_equations.dat.
"""

import matplotlib.pylab as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.externals import joblib
import pandas as pd
import subprocess
import re
import h5py
import os
import random

# import the cloud type file and set up the ranges
cl_fn = 'Cloud_type_data.csv'
cl_df = pd.DataFrame.from_csv(cl_fn) # pandas dataframe

# set up a default profile file and amend the content during each iteration
# of the values in the dataframe
no_sims = 100. # the number of simulations to do, will be spread out.
levels = 10 # the number of level intervals over which to spread the cloud
tot_cld = np.float(cl_df.sum()['perc'])
cl_df['sims'] = np.round(cl_df['perc']/tot_cld*no_sims)
wvc = 40.
aot = 0.2
press = 950. #hPa
sza = 40.
phi0 = 90.
umu = 1.
phi = 0.
st_wl = 540.
en_wl = 755.

# number of PC's to keep and plot
n_comp = 5

# no of random cloud spectrums to plot
n_clds = 5

# libradtran input files
inp_default = 'data_files_path /usr/local/share/libRadtran/data/\n\
    output_quantity reflectivity \nmol_abs_param reptran\n\
    atmosphere_file tropics \nrte_solver disort \naerosol_default \
    \naerosol_species_file continental_average \
    \naerosol_set_tau_at_wvl 550 %.3f \nmol_modify H2O %.3f MM \
    \nwavelength 400.0 1400.0 \npressure %.3f\
    \nsource solar /usr/local/share/libRadtran/data/solar_flux/atlas_plus_modtran\
    \nsza %.6f \nphi0 %.6f \numu %.6f \nphi %.6f \
    \nalbedo 0.0 \naltitude 0.0 \nzout TOA \nwavelength %.8f %.8f \
    \noutput_user lambda uu \n%s_file 1D \
    %s \ncloudcover %s 1'

temp_df = pd.DataFrame()
# the hierarchical index to use
col_index1 = [] # Abbreviation
col_index2 = [] # Cloud top
ser_list = []

cl_prof_fn = 'cloud_temp.dat' # file name of profile file
for r in cl_df.iterrows(): # loop over cloud types
    cl_abbr = r[1].name
    cl_bot = r[1]['alt_bot']
    cl_top_min = r[1]['alt_top_min']
    cl_top_max = r[1]['alt_top_max']
    sims = r[1]['sims']
    wi = 'wc' if r[1]['iw'] == 'W' else 'ic'
    cl_top_range = np.linspace(cl_top_min, cl_top_max, sims, endpoint=True)
    for h in cl_top_range: # loop over ranges of cloud tops heights
        print 'Simulating %s cloud type at altitude %.3f m.' % (cl_abbr, h)
        hs = np.linspace(h, cl_bot, levels, endpoint=True) / 1000 # convert to km's
        lwc = np.repeat(r[1]['lwc'], levels)
        rad = np.repeat(r[1]['rad'], levels)
        cl_arr = np.array([hs, lwc, rad]).T
        np.savetxt(cl_prof_fn, cl_arr, delimiter=' ', fmt='%.5f')
        inp_cloud = inp_default \
        %(aot, wvc, press, sza, phi0, umu, phi, st_wl, en_wl+1, wi, \
            cl_prof_fn, wi)
        process = subprocess.Popen('uvspec', stdin=subprocess.PIPE, stdout=\
            subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        spectrum, err = process.communicate(input=inp_cloud)
        spectrum = re.split('[\n\s]+', spectrum)[1:-1]
        spectrum = np.array(map(float, spectrum))
        spectrum = np.reshape(spectrum, (-1,2))
        lam = spectrum[:,0]
        spectrum = spectrum[:,1]
        spec_ser = pd.Series(data=spectrum, index=lam)
        ser_list.append(spec_ser)
        col_index1.append(cl_abbr)
        col_index2.append(h)

# save to Dataframe and then hdf5 file on disk
spec_df = pd.DataFrame(data=ser_list).transpose()
mul_index = pd.MultiIndex.from_arrays([col_index1, col_index2], \
    names=['type','cl_top_max'])
spec_df.columns = mul_index
spec_df.index.name = 'wavelength'

h5_fn = 'cloud_sims.h5'
if os.path.isfile(h5_fn):
    os.remove(h5_fn)
store = pd.HDFStore(h5_fn)
store['reflectance'] = spec_df

# wavelenghts for plotting
wls = spec_df.index

# plot the cloud spectrum, uncomment to use

rand_list = random.sample(np.arange(len(spec_df.columns)), n_clds)
cld_sub = spec_df[rand_list]
cld_sub.plot(title='Cloud Reflectance')
plt.ylabel('Reflectance')
plt.ylim([0,1.1])
plt.xlim(min(wls), max(wls))

# fit the PCA to the cloud data
cloud_specs = np.array(spec_df).T
pca = PCA(n_components=n_comp)
pca.fit(cloud_specs)

# plot the PC's
# PC's could be negative if the cloud reflectances reduce in reflectance from
# 1st to last, and vice versa for positives.
fig, axes = plt.subplots(nrows=pca.n_components_, ncols=1, sharex=True)
fig.set_figheight(10)
fig.set_figwidth(8)
plt.suptitle('Principle Components', y=1)
for i, (comps, ax, var) in enumerate(zip(pca.components_, axes, \
    pca.explained_variance_ratio_)):
    ax.plot(wls, comps)
    ax.set_ylabel('PC # %d\nvar. %8.5f%%' %(i+1, var*100))
    ax.set_xlim(min(wls), max(wls))
plt.xlabel('Wavelengths (nm)')
plt.tight_layout()
plt.show()

# save the pca into a dataframe and onto disk file
pc_cols1 = ['pc'+str(i) for i in range(1, n_comp+1)]
pc_cols2 = [i for i in pca.explained_variance_ratio_]
pca_df = pd.DataFrame(pca.components_.T, index=wls, columns=[pc_cols1, \
    pc_cols2])
pca_df.index.name = 'wavelength'
pca_df.columns.names = ['PCs', 'explained_variance_ratio']

store['pca'] = pca_df

# an alternative method to plot the PCs using pandas
# uncomment below and comment out the plotting method above
'''
axes = pca_df.plot(subplots=True, sharex=True, title='Principle Components', \
    legend=False)
for ax, label1, label2 in zip(axes.ravel(), pca_df.columns.levels[0], \
        pca_df.columns.levels[1]):
    label = '%s\n%8.5f%%' % (label1, label2*100)
    ax.set_ylabel(label)
    ax.set_xlabel('')
plt.xlabel('wavelength (nm)')
#plt.ylabel(pca_df.columns.levels[0])
'''
# save the pca object to a pickle file
pkl_fn = 'pca_object.p'
joblib.dump(pca, pkl_fn)

# print the information contents
print 'Explained variance in % per PC:'
for i, p in enumerate(pca.explained_variance_ratio_):
    print 'PC # %d: %8.5f%%' %(i+1, p*100)

store.close()