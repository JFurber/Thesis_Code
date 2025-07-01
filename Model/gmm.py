#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 13:12:33 2024

@author: jf01028
"""

"""
Useful Resources:
https://github.com/bnsreenu/python_for_microscopists/blob/master/52b-Understanding-GMM.py
https://github.com/scikit-learn/scikit-learn/blob/main/examples/mixture/plot_gmm_pdf.py
https://github.com/brianspiering/gaussian_mixture_models
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
http://www.cs.man.ac.uk/~fumie/tmp/bishop.pdf
https://scikit-learn.org/stable/modules/mixture.html
https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall06/reading/mixtureModels.pdf
http://ethen8181.github.io/machine-learning/clustering/GMM/GMM.html
https://www.cs.toronto.edu/~urtasun/courses/CSC411_Fall16/13_mog.pdf
https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
https://youtu.be/6gUdlygtscI
https://brilliant.org/wiki/gaussian-mixture-model/#:~:text=A%20Gaussian%20mixture%20model%20is,component%20means%20and%20variances%2Fcovariances.

"""

"""
Requires data to be shape (x,2), where x is the amount of data points

"""

from sklearn.mixture import GaussianMixture as GMM
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.patches import Ellipse
import seaborn as sns; sns.set()
from math import sqrt, log, exp, pi
from random import uniform
import pandas as pd
from scipy.io import savemat


def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    ax = ax or plt.gca()
    
    # Convert covariance to principal axes
    if covariance.shape == (2, 2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    
    # Draw the Ellipse
    for nsig in range(1, 4):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))
        
def plot_gmm(gmm, X, label=True, ax=None):
    ax = ax or plt.gca()
    labels = gmm.fit(X).predict(X)
    if label:
        ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', zorder=2)
    else:
        ax.scatter(X[:, 0], X[:, 1], s=40, zorder=2)
    ax.axis('equal')
    
    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha=w * w_factor)
        
def plot_gmm1(gmm, X, label=True, ax=None):
    ax = ax or plt.gca()
    # labels = gmm.fit(X).predict(X)
    if label:
        ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')
    else:
        ax.scatter(X[:, 0], X[:, 1], s=40)
    ax.axis('equal')
    
    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
        draw_ellipse(pos, covar, alpha=w * w_factor)
        
        
       
        

#%%

# data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data From Dez/Data/Woodchester_UTM.csv')
# data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data from Cornwall/Data/Cornwall_C2.csv')
# data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data from Cornwall/Data/C4_filtered_UTM.csv')
# data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data from Cornwall/Data/F1_filtered_UTM.csv')
# data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data from Cornwall/Data/Cornwall_F2.csv')
data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/0.Northern Ireland/Data/Ireland_AUGUST2024.csv')

x = data['Easting0']
x = pd.Series.to_numpy(x)
# x = x - np.min(x)

y = data['Northing0']
y = pd.Series.to_numpy(y, dtype='float64')
# y = y - np.min(y)

#Location data:
X = np.vstack((x,y)) 

del x, y, data

X_data = X.T

plt.scatter(*X)

xmin = np.min(X[0,:])
ymin = np.min(X[1,:])
xmax = np.max(X[0,:])
ymax = np.max(X[1,:])

#%%

# %matplotlib inline
# %matplotlib qt


# np.random.seed(4) #Cornwall, Woodchester
# np.random.seed(4) #NI

num_covar = 11

gmm = GMM(n_components= num_covar).fit(X_data)
labels = gmm.predict(X_data)

#generate the parameters
means = gmm.means_
covariances = gmm.covariances_
weights = gmm.weights_


# Store these matrices in a dictionary



# Save the dictionary to a .mat file
savemat('covariances.mat', {f'matrix_{i+1}': covariances[i] for i in range(num_covar)})

savemat('means.mat', {'means': means})

savemat('weights.mat',{'weights': weights})


#plot the data with the colour allocation of the clusters
plt.scatter(X_data[:, 0], X_data[:, 1], c=labels, s=40, cmap='viridis');        

# plot the data with the elipses underneath
plot_gmm(gmm, X_data)

plot_gmm1(gmm, X_data)


#%% Negative log-likelihood predicted by a GMM (i.e. contours)

# g_both = [gmm.pdf(e) for e in x]
x = np.linspace(-0.1,4500,1000)
y = np.linspace(-0.1, 2750,1000)
X, Y = np.meshgrid(x, y)
XX = np.array([X.ravel(), Y.ravel()]).T

Z = -gmm.score_samples(XX)
Z = Z.reshape(X.shape)

CS = plt.contour(
    X, Y, Z, norm=LogNorm(vmin=1.0, vmax=100.0), levels=np.logspace(0, 1, 10)
)
CB = plt.colorbar(CS, shrink=0.8, extend="both")
plt.scatter(X_data[:, 0], X_data[:, 1], 0.8)
plt.title("Negative log-likelihood predicted by a GMM")
plt.axis("tight")
plt.show()

#%% 3d plane of above log-likelihood prediction

fig = plt.figure(figsize=(13, 7))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='coolwarm', edgecolor='none',vmin=13, vmax=25)
ax.set_xlabel('x')
ax.set_ylabel('y') 
ax.set_zlabel('PDF')
# ax.set_zlim([0, 15])
ax.set_title('Surface plot of Gaussian 2D KDE')
fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF
# ax.view_init(60, 35)


