#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 08:07:21 2024

@author: jf01028
"""


import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt

# Define the Gaussian Mixture Model
class GMM:
    def __init__(self, means, covariances, weights):
        self.means = means
        self.covariances = covariances
        self.weights = weights

    def potential(self, x):
        return -np.log(sum(w * multivariate_normal(mean=m, cov=c).pdf(x) for m, c, w in zip(self.means, self.covariances, self.weights)))



class SDE:
    def __init__(self, dim, dt, gmm, interaction_strength,S):
        
        self.dim = dim
        self.dt = dt
        self.gmm = gmm
        self.interaction_strength = interaction_strength
        self.eps = 1e-5  # Small value to compute numerical gradient
        self.S = S
        self.n = S.shape[0]
    
    # def Attraction(self, x, i, j): 
    #     return -2*(x[2*i:2*i+2] - x[2*j:2*j+2])
    
    def Attraction(self, x, i, j): 
        xi = x[2*i:2*i+2]
        xj = x[2*j:2*j+2]
        # Debugging print statements
        # print(f"i: {i}, j: {j}, xi shape: {xi.shape}, xj shape: {xj.shape}")
        if xi.shape[0] != 2 or xj.shape[0] != 2:
            raise ValueError(f"Incorrect slicing for agents i={i}, j={j}. Shapes are xi: {xi.shape}, xj: {xj.shape}")
        return -2*(xi - xj)
    
    # def Repulsion(self, x, i, j):
    #     return np.array([ 2*(x[2*i] - x[2*j])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2, 2*(x[2*i+1] - x[2*j+1])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2])

    def Repulsion(self, x, i, j, epsilon=1e-8):
        dx = x[2*i] - x[2*j]
        dy = x[2*i+1] - x[2*j+1]
        denominator = (dx**2 + dy**2 + epsilon)**2  # Add epsilon to avoid division by zero (this would occur when two agents are at the same position)
        return np.array([2*dx/denominator, 2*dy/denominator])
    
    def drift(self,x):
        y = np.zeros_like(x,dtype='float64') # Returns an array of given shape and type as the given array (i.e. x), with zeros
        for i in range(0, self.n):
            position_i = x[2*i:2*i+2]
            gradV = np.zeros(self.dim)
            for k in range(self.dim):
                dx = np.zeros(self.dim)
                dx[k] = self.eps
                gradV[k] = (self.gmm.potential(position_i + dx) - self.gmm.potential(position_i - dx)) / (2 * self.eps)
            # dx = np.zeros(self.dim)
            # dx[i] = self.eps
            # gradV = (self.gmm.potential(position_i + dx) - self.gmm.potential(position_i - dx)) / (2 * self.eps)
            y[2*i:2*i+2] = -gradV
            for j in range(0, self.n):
                if self.S[i][j] == 1:
                    attraction = self.interaction_strength * self.Attraction(x, i, j)
                    y[2*i:2*i+2] += attraction
                elif self.S[i][j] == -1:
                    repulsion = self.interaction_strength * self.Repulsion(x, i, j)
                    y[2*i:2*i+2] += repulsion
        return y
    
    def diffusion(self):
        return 22* np.sqrt( self.dt)


    def step(self, x):
        # print(f"Before step, x: {x[:4]}")  # Print first few positions for brevity
        drift = self.drift(x)                                            
        diffusion = self.diffusion() * np.random.normal(size=x.shape[0])
        x_new = x + drift * self.dt + diffusion
        # print(f"After step, x: {x_new[:4]}")
        return x_new

num_agents = 10  # Number of agents
num_steps = 100000  # Number of steps in the simulation

# S = np.ones((2,2))
# S = np.array([[0,1,1],[1,0,1],[1,1,0]])

S1=2*np.random.randint(0,2,size =(num_agents,num_agents))-1 #randomised attraction/repulsion matrix
np.fill_diagonal(S1, 0)
# print(S)

# Create the GMM
gmm = GMM(means, covariances, weights) #Import means, covariances and weights from data


steps = np.zeros((num_steps, 2 * num_agents))
# x = np.array([ 600,1500,600,1525,625,1517],dtype='float64')
# x = np.array([500, 1000,1000, 2000, 2000, 500, 3000, 2000,500, 1000,1000, 2000, 2000, 600, 3000, 2000, 600,1500,600,2000],dtype='float64') #
x = np.array([850,1420,900,1420,3870,580,3410,1710,3810,1350,1720,1300,2950,680,1600,2020,2690,1870,1500,1200],dtype='float64') # start near means


sde = SDE(2, 0.01, gmm,0.0001,S1)  # Each agent is 2D


import time

start = time.time()


for step in range(num_steps):
    x = sde.step(x)  # Update x with the new state returned by the step method
    steps[step,:]=x


        
end = time.time()
print(end - start)

plt.figure()
plt.plot(*steps.T)
plt.scatter(means[:,0],means[:,1],marker='x')




# steps = steps.T

# plt.figure()
# plt.clf()
# for i in range(num_agents):
#     plt.plot(steps[2*i,:], steps[2*i+1,:])
# plt.xlabel('Easting')
# plt.ylabel('Northing')


# plt.figure()
# plt.plot(steps[0,:], steps[1,:])



#%% Plot of Results - create the base
import d3s.domain as domain

xmin = np.min(X[0,:])
ymin = np.min(X[1,:])
xmax = np.max(X[0,:])
ymax = np.max(X[1,:])
bounds = np.array([[xmin, xmax],[ymin,ymax]])

x1 = np.linspace(xmin,xmax,100)
y1 = np.linspace(ymin,ymax,100)
x1, y1 = np.meshgrid(x1, y1)


#Set how many boxes you the potential to be visualised over. 
#The larger the number the longer it takes, the smaller the number, the less accurate the potential.
boxes = np.array([100, 100])

Omega = domain.discretization(bounds, boxes) 

PsiC = KDE.V(Omega.midpointGrid()) #potential
PsiC_KDE1 = np.reshape(PsiC, (100, 100))

#%% 



#Contour Plot
plt.figure()
plt.clf()
Omega.plot(PsiC[0,:], '2D') #colour plot of energy potential
# plt.contour(x1, y1, PsiC_KDE1,100)
# plt.clim(3000,6000) #colour limits for the colourplot. Check 3D plot for best limit to use.
for i in range(num_agents):
    plt.plot(steps_kde.T[2*i,:], steps_kde.T[2*i+1,:], color="black")
    plt.plot(steps_kde.T[2*i,0], steps_kde.T[2*i+1,0],'x')

plt.xlabel('Easting')
plt.ylabel('Northing')



#Save to file and plot in MATLAB (check folder directory first to where it is being saved)
# np.savetxt('14M_traj_01_10000.csv',X_new1.T,delimiter = ',')


plt.figure()
j=0
plt.plot(steps_kde1.T[j,:], steps_kde1.T[j+1,:], color="black")
j=1
plt.plot(steps_kde1.T[j,:], steps_kde1.T[j+1,:], color="black")


