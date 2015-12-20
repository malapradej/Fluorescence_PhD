'''A script in which a selection of a layers elements are made
prior to running this script. The script is then modified to 
choose from 3 different features (x, y and z)which are plotted 
by this script in a scatterplot. To select only those elements
you want use select by expression or another.
==================================
input: none. Only change the feature_x etc text below.
output: scatterplot of the features attributes.
'''
feature_x = 'p'
feature_y = 'p_easy'
feature_z = 'clf_prior'

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

layer = qgis.utils.iface.activeLayer()
selected_features = layer.selectedFeatures()

x = []
y = []
z = []

for f in selected_features:
    x.append(f[feature_x])
    y.append(f[feature_y])
    z.append(f[feature_z])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
#plt.scatter(x, y)
ax.set_xlabel(feature_x)
ax.set_ylabel(feature_y)
ax.set_zlabel(feature_z)
plt.show()

