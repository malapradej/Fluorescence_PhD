#!/usr/bin/python
"""
A script to do a PCA of all the LUTs in the different view, solar zenith and
relative azimuth angle directories. The script will loop over every folder and
create a PC object file in pickled format and PC's in hdf5 format with all the
detail which can then be loaded into our ipyton notebook (WB2) for further
analysis.
"""
import numpy as np
import pandas as pd
import cPickle
import os
from sklearn.externals import joblib
from sklearn.decomposition import PCA
import glob

# select the wavelengths range over the red-edge
min_wl = 680
max_wl = 755

# select the number of PCs to calculate
n_comp = 10

# loop over all directories starting in the defined root and do the pca
# on all the picled files encountered.
rootdir = '/home/malapradej/Documents/PhD_UCL/Data/LUT/LUT3'
for dirname, subdirs, files in os.walk(rootdir):
    if len(subdirs) > 0:
        continue
    print 'PCs of %s.' % dirname
    pattern = os.path.join(dirname, '*_*_*_*.p')
    files = glob.glob(pattern)
    for i, fn in enumerate(files):
        dic = cPickle.load(open(fn, 'rb'))

        # set up spectral dataframes for all files in folder at 1st file
        if i == 0:
            lam = dic['lam']
            ind = np.logical_and(lam>=min_wl, lam<=max_wl)
            lam = lam[ind]
            a_df = pd.DataFrame(index=lam)
            d_df = a_df.copy()
            s_df = a_df.copy()

        atm_path = dic['atm_path'][ind]
        dbl_trans = dic['dbl_trans'][ind]
        spher_alb = dic['spher_alb'][ind]

        a_df[fn] = atm_path
        d_df[fn] = dbl_trans
        s_df[fn] = spher_alb

    # different PCAs fitted and saved to disk
    pca_a = PCA(n_components=n_comp)
    pca_a.fit(np.array(a_df.T))
    fn = os.path.join(dirname, 'pca_atm_path.p')
    joblib.dump(pca_a, fn)

    pca_d = PCA(n_components=n_comp)
    pca_d.fit(np.array(d_df.T))
    fn = os.path.join(dirname, 'pca_dbl_trans.p')
    joblib.dump(pca_d, fn)

    pca_s = PCA(n_components=n_comp)
    pca_s.fit(np.array(s_df.T))
    fn = os.path.join(dirname, 'pca_spher_alb.p')
    joblib.dump(pca_s, fn)

    # create dataframe files of all PCs and explained var ratios and save to disk
    pc_cols1 = ['PC'+str(i).zfill(2) for i in range(1, n_comp+1)]

    pc_cols2 = [i for i in pca_a.explained_variance_ratio_]
    pca_a_df = pd.DataFrame(pca_a.components_.T, index=lam, columns= [pc_cols1, \
        pc_cols2])

    pc_cols2 = [i for i in pca_d.explained_variance_ratio_]
    pca_d_df = pd.DataFrame(pca_d.components_.T, index=lam, columns= [pc_cols1, \
        pc_cols2])

    pc_cols2 = [i for i in pca_s.explained_variance_ratio_]
    pca_s_df = pd.DataFrame(pca_s.components_.T, index=lam, columns= [pc_cols1, \
        pc_cols2])

    fn = os.path.join(dirname, 'pca.h5') # remove if file exists
    if os.path.exists(fn):
        os.remove(fn)
    store = pd.HDFStore(fn, 'a')
    store['atm_path_pca'] = pca_a_df
    store['dbl_trans_pca'] = pca_d_df
    store['spher_alb_pca'] = pca_s_df
    store.close()

    # plot the PCs
    '''
    pca_a_df.plot(subplots=True, sharex=True, figsize=(10,10), \
        title='PCA atm_path')
    pca_d_df.plot(subplots=True, sharex=True, figsize=(10,10), \
        title='PCA dbl_trans')
    pca_s_df.plot(subplots=True, sharex=True, figsize=(10,10), \
        title='PCA spher_alb')
    '''