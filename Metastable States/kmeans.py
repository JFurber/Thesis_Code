#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code is used to determine the amount of 'k' clusters, using k-means algorithm.

Firstly import the data and separate the long/lat variables.

Next, estimate the amount of clusters in the data, this is done using the elbow plot method and the silhouette plot method.

Once chosen, split the data into k amount of clusters.

Data can be exported along the way to generate figure outside of Python.
"""

from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from yellowbrick.cluster import KElbowVisualizer
from yellowbrick.cluster import SilhouetteVisualizer


# Test Data:

# x1 = np.array([3, 1, 1, 2, 1, 6, 6, 6, 5, 6, 7, 8, 9, 8, 9, 9, 8])
# x2 = np.array([5, 4, 5, 6, 5, 8, 6, 7, 6, 7, 1, 2, 1, 2, 3, 2, 3])
# X = np.vstack((x1,x2))
# X = X.T

# Go to cell 'Elbow Plot' from here

# %% Import data

#Import Data - full data set
data = pd.read_csv('filename.csv',low_memory=False) #paste path file here

# Filter for a specific year - this is not a requirement if using the full data set.
data = data[data['year']==2022]

# Pull out the coordinate data - ensure you have the correct column name
x = data['Easting']
x = pd.Series.to_numpy(x)
# x = x - np.min(x) #If you want to move the coordinates to a (0,0) base, uncomment this line

y = data['Northing']
y = pd.Series.to_numpy(y, dtype='float64')
# y = y - np.min(y) #If you want to move the coordinates to a (0,0) base, uncomment this line

#Location data:
X = np.vstack((x,y))
X = X.T

#Visualise Data
plt.plot()
plt.scatter(*X.T)
plt.show()

del x,y,data

#For ease and saving along the way, name site here:
Site = 'Site'

#Save the coordinates separately for plotting
np.savetxt(f'{Site}_coords.csv',X,delimiter=',')

#%% Elbow Plot

km = KMeans(random_state=42,n_init=100)
visualizer = KElbowVisualizer(km, k=(2,11)) #generate a distortion score for 2 clusters up to 10 clusters

visualizer.fit(X) # Fit the data to the visualizer
visualizer.show() # Plots elbow plot, with fit time in seconds. Plots where the algorithm believes the elbow is at.

# Save variables to plot out side of Python (if required)
np.savetxt(f'{Site}_k_scores.csv',visualizer.k_scores_)
np.savetxt(f'{Site}_k_times.csv',visualizer.k_timers_)
np.savetxt(f'{Site}_k_vaules.csv',visualizer.k_values_)

#%% Silhouette Plots

k = 3 # number of predicted clusters

#Print 1 silhouettes
model = KMeans(k, random_state=42, n_init=100)
visualizer_1SP = SilhouetteVisualizer(model, colors='yellowbrick')

visualizer_1SP.fit(X)        # Fit the data to the visualizer
visualizer_1SP.show()        # Finalize and render the figure

#print multiple silhouettes - uncomment lines below if data is to be saved outside of python
for i in [2,3, 4]: # Can choose suitable range chosen based on elbow score
    model = KMeans(n_clusters = i, n_init=100)
    cluster_labels = model.fit_predict(X)
    visualizer_SP = SilhouetteVisualizer(model)
    # np.savetxt(f'{Site}_silhouette_labels_{i}.csv',cluster_labels )
    visualizer_SP.fit(X)    # Fit the data to the visualizer
    y = visualizer_SP.silhouette_score_
    # f = open(f'{Site}_sil_avg_{i}.txt', "w") 
    # f.write(str(y))
    # f.close()
    visualizer_SP.show()

#%% Generate clusters

k = 3 # Chosen amount of clusters

kmeans = KMeans(n_clusters = k,init='k-means++',n_init=300, max_iter = 500)

# Predict the labels of clusters and save
label = kmeans.fit_predict(X)
np.savetxt(f'{Site}_label_k{k}.csv',label,delimiter=',')

#Getting the Centroids
centroids = kmeans.cluster_centers_
np.savetxt(f'{Site}_centroids_k{k}.csv',centroids,delimiter=',')

#plotting the results (clusters only):
u_labels = np.unique(label)
for i in u_labels:
    plt.scatter(X[label == i , 0] , X[label == i , 1] , label = i)
plt.scatter(centroids[:,0] , centroids[:,1] , s = 80, color = 'k')
plt.legend()
plt.show()


