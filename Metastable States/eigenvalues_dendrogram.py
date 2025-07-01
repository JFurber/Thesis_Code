#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 17:11:20 2024

@author: jf01028
"""

# clustering dataset
# determine k using elbow method

import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl #needed for VSCode


#Import Data - full data set

X = pd.read_csv('.csv',header=None) #paste path file here


X = X.to_numpy()


#%

import numpy as np
from sklearn.cluster import KMeans
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

data = X.reshape(-1)


# Optimal number of clusters using dendrogram (cut at appropriate height)
from scipy.cluster.hierarchy import fcluster
max_d = 0.001  # Maximum distance for cluster separation
hierarchical_labels = fcluster(sch.linkage(data.reshape(-1, 1), method='ward'), max_d, criterion='distance')

print("Hierarchical Clustering Labels:", hierarchical_labels)

n = len(X)

plt.scatter(np.linspace(0,n-1,n),X)
plt.show()


plt.scatter(np.linspace(0,n-1,n),X,c=hierarchical_labels,cmap='ocean')
plt.show()


linkage_matrix = sch.linkage(data.reshape(-1, 1), method='ward')
sch.dendrogram(linkage_matrix)
# Hierarchical Clustering
dendrogram = sch.dendrogram(sch.linkage(data.reshape(-1, 1), method='ward'))
plt.title('Dendrogram')
plt.xlabel('Data Points')
plt.ylabel('Euclidean Distances')
plt.show()

leaf_order = dendrogram['leaves']

np.savetxt("linkage_matrix.csv", linkage_matrix, delimiter=",")
np.savetxt("Hierarchial_labels.csv", hierarchical_labels, delimiter=",")
np.savetxt("eigenvalues.csv",data,delimiter=",")
np.savetxt("leaf_order.csv", leaf_order, delimiter=",")

#%%

# x1 = np.array([3, 1, 1, 2, 1, 6, 6, 6, 5, 6, 7, 8, 9, 8, 9, 9, 8])
# x2 = np.array([5, 4, 5, 6, 5, 8, 6, 7, 6, 7, 1, 2, 1, 2, 3, 2, 3])

# plt.plot()
# plt.xlim([0, 10])
# plt.ylim([0, 10])
# plt.title('Dataset')
# plt.scatter(x1, x2)
# plt.show()

x = np.linspace(0, 19,20)

# plt.scatter(x,X.T)

# create new plot and data
# plt.plot()
# X = np.array(list(zip(x1, x2))).reshape(len(x1), 2)
# colors = ['b', 'g', 'r']
# markers = ['o', 'v', 's']

# k means determine k
distortions = []
K = range(1,10)
for k in K:
    kmeanModel = KMeans(n_clusters=k).fit(X)
    kmeanModel.fit(X)
    distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclid'), axis=1)) / X.shape[0])
# np.savetxt('kmeans_woodchester1.csv',distortions,delimiter=',')

# Plot the elbow
plt.plot(K, distortions, 'bx-')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
# plt.savefig('elbow_woodchester1.jpg')




