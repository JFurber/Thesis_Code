#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Within this code we want to estimate the energy potential from data using Kernel Density Estimation (KDE).

1. Import data. Require the paired co-ordinates (preferably in UTM format). Translate the coordinates to (0,0) to keep anonymity of the data

2. Next, it is ideal to estimate a value for the bandwidth h. This value shall be used within KDE as the smoothing parameter. 
The link below gives the details of three estimators that can be used (normal reference, cross-validation least squares, and cross-validation maximum likelihood).

https://www.statsmodels.org/devel/generated/statsmodels.nonparametric.kernel_density.KDEMultivariate.html

3. Generate the potential using KDE from the toolbox created from Stefan Klus (found on Github). Save the spydata for later use.

https://github.com/sklus/d3s

4. Visualise the potential - either on Matlab or Python

"""

import numpy as np
import pandas as pd
import os
import statsmodels.api as sm
import matplotlib.pyplot as plt
import d3s.domain as domain
import d3s.kernels as kernals

import time, sys

#%% Step 1: Import the Data

start = time.time()

site = 'C4'

#Import Data - full data set
data = pd.read_csv('file.csv') #paste path file here

x = data['Easting']
x = pd.Series.to_numpy(x)
x = x - np.min(x)

y = data['Northing']
y = pd.Series.to_numpy(y, dtype='float64')
y = y - np.min(y)

#Location data:
X = np.vstack((x,y)) 

del x,y, data #clean variable explorer

print(f'Cornwall {site}')
sys.stdout.flush()

#%% Step 2: Estimate a value for the bandwidth
#Option of three estimators. The lower interval of normal reference can be a good choice is unable to choose.

# Normal Reference
# h_NR = sm.nonparametric.KDEMultivariate(data=X, var_type='cc', bw='normal_reference')
# print('Method NR gives:')
# print(h_NR.bw)
# sys.stdout.flush()

# # Cross-validation Least Squares
# h_cvls = sm.nonparametric.KDEMultivariate(data=X, var_type='cc', bw='cv_ls')
# print('Method CVLS gives:')
# print(h_cvls.bw)
# sys.stdout.flush()

# # Cross-validation Maximum Liklihood
# h_cvml = sm.nonparametric.KDEMultivariate(data=X, var_type='cc', bw='cv_ml')
# print('Method CVML gives:')
# print(h_cvml.bw)
# sys.stdout.flush()


# end = time.time()
# print(end - start)
#%% Step 3: Generate the potential
h_value = 143
noise = 25.015
h = kernals.gaussianKernel(143.0169388) # Chosen bandwidth, potentially estimated in step 2. 
D = (kernals.densityEstimate(X, h , beta = 2/(noise**2))) #KDE, with X the imported data and h the bandwidth. beta is the inverse temperature (for MD applications) so can be set to 1.


print(f' h {h_value}')
print(f' noise {noise}')

#%% Step 4: Visualise the potential

#Set bounds for the plot, i.e. minimum coordinates to the max in both directions with a gap
bounds = np.array([[-150, 4450],[-150,4650]])

print(bounds)

#Set how many boxes you the potential to be visualised over. 
#The larger the number the longer it takes, the smaller the number, the less accurate the potential.
boxes = np.array([250, 250])

Omega = domain.discretization(bounds, boxes) 

PsiC = D.V(Omega.midpointGrid()) #potential
print('PsiC computed')
sys.stdout.flush()

C = Omega.midpointGrid()
print('C computed')
sys.stdout.flush()

# plt.figure()
# Omega.plot(PsiC[0,:], '3D')  #plot the potential in Python


PsiC1 = D.rho((Omega.midpointGrid())) #invariant density
print('PsiC1 computed')
sys.stdout.flush()
# plt.figure(3)
# Omega.plot(PsiC1[0,:], '3D')

#Save the data: (check folder directory first to where it is being saved) to plot in Matlab
np.savetxt(f'potential_{site}_h{h_value}.csv',PsiC.T,delimiter=',')
np.savetxt(f'density_{site}_h{h_value}.csv',PsiC1.T,delimiter=',')
np.savetxt(f'C_{site}.csv',C.T,delimiter=',')
 
end = time.time()
print(end - start)
