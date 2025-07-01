#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is used to determine the metastable states within a dynamical system.

It makes use of the code d3s, created by Stefan Klus, found at https://github.com/sklus/d3s

See Matlab file for interpolation of data. Generate 3 files:
    1. a file with the first to second to last coordinates (X)
    2. a file with the second to last coordinates (Y)
    3. a file with all the interpolated values (not split)

Then, import the data and separate the easting/northing variables.

Choose your bounds and boxes - here chosen a buffer with 150m, and boxes chosen to obtain 100m x 100m

Choose how many eigenvalues you want to compute. Completing k-means on the raw data (computing elbow plot/silhoutte plot can help guide this choice)

Observing the eigenvalue plot - identify the spectral gap:
    cluster the number of eigenfunctions based on this choice.
    
If you also want to calculate transition probabilities based on the metastable states, this is done at the end.

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import d3s.domain as domain
import d3s.observables as observables
import d3s.algorithms as algorithms

from sklearn.cluster import KMeans



#%% Import Data

dataX = pd.read_csv('X.csv', header=None) #paste path file here (file with first coordinates to second to last coordinates)
dataY = pd.read_csv('Y.csv',header=None) #paste path file here (file with second coordinates to last coordinates)

#For ease and saving along the way, name site here:
Site = 'Site'

x = dataX[0]
x = pd.Series.to_numpy(x)

y = dataX[1]
y = pd.Series.to_numpy(y, dtype='float64')

X = np.vstack((x,y))

x = dataY[0]
x = pd.Series.to_numpy(x)

y = dataY[1]
y = pd.Series.to_numpy(y, dtype='float64')

Y = np.vstack((x,y))

del x,y, dataY, dataX #clean variable explorer


#%% Set up EDMD

"""
Set the bounds for domain - here added a 150m buffer around the edge and rounded.
"""
bounds = np.array([[round(np.min(X[0,:])-150,-1),round(np.max(X[0,:])+150,-1)],[round(np.min(X[1,:])-150,-1), round(np.max(X[1,:])+150,-1)]])

"""
Choose the amount of boxes the domain will be split into. Can choose different values based on the domain (i.e. to obtain 100mx100m boxes)
"""
box_x = np.round((bounds[0,1]-bounds[0,0])/100).astype(int)
box_y = np.round((bounds[1,1]-bounds[1,0])/100).astype(int)
boxes = np.array([box_x, box_y])

Omega = domain.discretization(bounds, boxes)

C = Omega.midpointGrid()

# np.savetxt('C.csv',C, delimiter=',') # used for plotting outside of python

# Choose observables:
# psi = observables.monomials(10) #found monomials were not best for the GPS data
psi = observables.indicators(Omega)

"""
Choose number of eigenvalues you want to compute - the larger the number the longer it will take
"""
evs = 40 # number of eigenvalues/eigenfunctions to be computed

PsiC = psi(Omega.midpointGrid()) # observables evaluated at midpoints of the grid

#%% Run EDMD

_, d, V = algorithms.edmd(X, Y, psi, operator='K', evs=evs)

"""
Plot the eigenvalues
"""
plt.figure()
plt.plot(np.real(d),'.')
plt.xlabel('Eigenvalue Index')
plt.ylabel('Eigenvalue')
plt.show()

"""
Plot the first n eigenfunctions - expect the first one to be constant
"""
for i in range(1):
    plt.figure()
    r = np.real(V[:,i].T @ PsiC)
    Omega.plot(r, '2D')
    plt.title('EDMD K, eigenfunction %d' % i)
    plt.show()

np.savetxt(f'{Site}_eigenvalues15.csv',d, delimiter=',')
np.savetxt(f'{Site}_eigenvalues15_real.csv', np.real(d),delimiter=',')

#%% Find Spectral Gap

#Dendogram

import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster

eigenvalues = np.real(d).reshape(-1)

# Optimal number of clusters using dendrogram (cut at appropriate height)
max_d = 0.005  # Maximum distance for cluster separation
hierarchical_labels = fcluster(sch.linkage(eigenvalues.reshape(-1, 1), method='ward'), max_d, criterion='distance')

linkage_matrix = sch.linkage(eigenvalues.reshape(-1, 1), method='ward')
sch.dendrogram(linkage_matrix)

# Hierarchical Clustering
dendrogram = sch.dendrogram(sch.linkage(eigenvalues.reshape(-1, 1), method='ward'))
plt.title('Dendrogram')
plt.xlabel('Data Points')
plt.ylabel('Euclidean Distances')
plt.show()

n = len(eigenvalues)
plt.scatter(np.linspace(0,n-1,n),eigenvalues,c=hierarchical_labels,cmap='ocean')
plt.show()

np.savetxt('linkage_matrix.csv', linkage_matrix, delimiter=",")
np.savetxt('Hierarchial_labels.csv', hierarchical_labels, delimiter=",")

print("Hierarchical Clustering Labels:", hierarchical_labels)
unique, counts = np.unique(hierarchical_labels, return_counts=True)
dict(zip(unique, counts))

EV_diff = np.zeros(((n-1),1))
for i in range(n-1):
    EV_diff[i] = abs(eigenvalues[i]-eigenvalues[i+1])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(eigenvalues[0:41], 'g-')
ax2.plot(EV_diff[0:41], 'b-')

ax1.set_xlabel('X data')
ax1.set_ylabel('Y1 data', color='g')
ax2.set_ylabel('Y2 data', color='b')

plt.xticks([1, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41])  # Set x-axis gridline positions
ax1.grid(which='major', axis='both', linestyle='--')
ax2.grid()

plt.show()

np.savetxt('EV_diff.csv', EV_diff, delimiter=",")

plt.scatter(np.linspace(0,30,31),eigenvalues[0:31])

#%%
"""
Locate which boxes the coordinates are in, the kmeans clusters the boxes (not the coordinates).
"""

psiX = psi(X)

location_boxes = np.zeros((len(psiX[1])))
for i in range(len(psiX[1])):
    location_boxes[i] = np.where(psiX[:,i] == 1)[0][0]

"""
Find where there are no badgers in a box
"""

sum_psiX = np.sum(psiX,1)

index, = np.where(sum_psiX == 0)

#%%

"""
Organise eigenfunctions to complete kmeans

Only to need to compute the ones needed for kmeans, i.e. if you want to cluster for 4, only to need to compute r0-r3.

Comment out the ones not in use.
"""

a = np.arange(0,box_x*box_y)

r0 = np.real(V[:,0].T @ PsiC) #fist eigenfunction
r1 = np.real(V[:,1].T @ PsiC) #second eigenfunction
r2 = np.real(V[:,2].T @ PsiC) #third eigenfunction
r3 = np.real(V[:,3].T @ PsiC) #fourth eigenfunction
r4 = np.real(V[:,4].T @ PsiC) #fifth eigenfunction
r5 = np.real(V[:,5].T @ PsiC) #sixth eigenfunction
r6 = np.real(V[:,6].T @ PsiC) #seventh eigenfunction
r7 = np.real(V[:,7].T @ PsiC) #eighth eigenfunction
r8 = np.real(V[:,8].T @ PsiC) #ninth eigenfunction
r9 = np.real(V[:,9].T @ PsiC) #tenth eigenfunction
r10 = np.real(V[:,10].T @ PsiC) #eleventh eigenfunction
r11 = np.real(V[:,11].T @ PsiC) #twelth eigenfunction


# np.savetxt('r0.csv',r0, delimiter=',')

"""
remove where there are no badgers - these are the eigenfunctions that will be clustered
"""
a_new = np.delete(a,index)
r0_new = np.delete(r0,index)
r1_new = np.delete(r1,index)
r2_new = np.delete(r2,index)
r3_new = np.delete(r3,index)
r4_new = np.delete(r4,index)
r5_new = np.delete(r5,index)
r6_new = np.delete(r6,index)
r7_new = np.delete(r7,index)
r8_new = np.delete(r8,index)
r9_new = np.delete(r9,index)
r10_new = np.delete(r10,index)
r11_new = np.delete(r11,index)


del r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 

#%% k-means

"""
Stack the first k eigenfunctions that you want to cluster, based on the spectral gap found in the eigenvalues

"""

R_Kmeans = np.vstack((r0_new,r1_new,r2_new,r3_new,r4_new))
# , r5_new, r6_new ,r7_new , r8_new, r9_new, r10_new, r11_new

k = 5

kmeans = KMeans(n_clusters=k).fit(np.real(R_Kmeans.T)).fetch_model()
c = kmeans.transform(R_Kmeans.T)


"""
Allocate the boxes/clusters to the coordinates in X
"""

c_coordinates = np.zeros((len(location_boxes)))
for i in range(len(location_boxes)):
    y = location_boxes[i]
    x = np.where(a_new == y)[0][0]
    c_coordinates[i] = c[x]

np.savetxt(f'{Site}_c_coordinates{k}.csv',c_coordinates) 

"""
Plot the coordinates with the clusters
"""
plt.scatter(*X, c = c_coordinates)
plt.show()

#%%
"""
To put interpolated data into the clusters in order to create Markov Chain

In other words, this allocates each coordinate with a metastable state, to go on to calculate transition rates between clusters.
"""

dataX = pd.read_csv('interpolate.csv', header=None) #paste path file here - these are the full interpolated coordinates (not split into X and Y files)

x = dataX[2]
x = pd.Series.to_numpy(x,dtype = 'float64')

y = dataX[3]
y = pd.Series.to_numpy(y, dtype='float64')

X_interp = np.vstack((x,y))

id = dataX[0]
id = pd.Series.to_numpy(id, dtype='int')

interval = dataX[1]
interval = pd.Series.to_numpy(interval)

dataX_df = np.vstack((id,interval,X_interp))

del x,y #clean variable explorer

psiX_interp = psi(X_interp)

location_boxes_interp = np.zeros((len(psiX_interp[1])))
for i in range(len(psiX_interp[1])):
    location_boxes_interp[i] = np.where(psiX_interp[:,i] == 1)[0][0]

"""
Find where there are no badgers in a box
"""

sum_psiX_location = np.sum(psiX_interp,1)
index_interp, = np.where(sum_psiX_location == 0)

a = np.arange(0,box_x*box_y) #625 boxes
a_new_interp = np.delete(a,index_interp)

c_coordinates_interp = np.zeros((len(location_boxes_interp)))
for i in range(len(location_boxes_interp)):
    y = location_boxes_interp[i]
    x = np.where(a_new_interp == y)[0][0]
    c_coordinates_interp[i] = c[x]

dataX_location = np.vstack((dataX_df, location_boxes_interp, c_coordinates_interp))
np.savetxt(f'{Site}_data_location.csv', dataX_location.T,delimiter = ',') 


