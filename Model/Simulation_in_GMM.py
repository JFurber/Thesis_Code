#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:15:08 2024

@author: jf01028
"""

# This also works!

import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt


# from sklearn.mixture import GaussianMixture

# Assume data is a 2D numpy array where each row is a sample
# data = np.random.rand(1000, 2)  # Replace this with your actual data

#To plot the contours:
def Gaussian_Mixture(x, y):
    return np.sum([weights[i] * multivariate_normal.pdf(np.dstack([x, y]), mean=means[i], cov=covariances[i]) for i in range(len(weights))], axis=0)


class GMM:
    def __init__(self, means, covariances, weights):
        self.means = means
        self.covariances = covariances
        self.weights = weights

    def potential(self, x):
        return -np.log(sum(w * multivariate_normal(mean=m, cov=c).pdf(x) for m, c, w in zip(self.means, self.covariances, self.weights)))
    
    
# Define the SDE
class SDE:
    def __init__(self, dim, dt, gmm):
        self.dim = dim #dimension
        self.dt = dt #time step
        self.gmm = gmm
        self.eps = np.sqrt(np.finfo(float).eps)

    def multivariate_gaussian_pdf(self, xv):
        return multivariate_normal.pdf(x, means, covariances)

    # Responsibilities (posterior probabilities) w_k(x)
    def compute_responsibilities(self, x):
        K = len(weights)
        responsibilities = np.zeros(K)
        total = sum(weights[k] * self.multivariate_gaussian_pdf(x, means[k], covariances[k]) for k in range(K))
        for k in range(K):
            responsibilities[k] = weights[k] * self.multivariate_gaussian_pdf(x, means[k], covariances[k]) / total
        return responsibilities

    def gradient_E(self, x): #works when only 1 agent
        K = len(weights)
        D = len(x)
        grad = np.zeros(D)
        responsibilities = self.compute_responsibilities(x, weights, means, covariances)
        for k in range(K):
            diff = x - means[k]
            inv_cov_k = np.linalg.inv(covariances[k])
            grad += responsibilities[k] * np.dot(inv_cov_k, diff)
        return -grad  # Gradient should be negative for minimization

    def drift(self, x): #Computes drift using numerical differentiation
        grad = np.zeros(self.dim)
        for i in range(self.dim):
            dx = np.zeros(self.dim)
            dx[i] = self.eps
            grad[i] = (self.gmm.potential(x + dx) - self.gmm.potential(x - dx)) / (2 * self.eps)
        return -grad

    def diffusion(self):
        return 22 * np.sqrt(self.dt)

    def step(self, x):
        return x + self.drift(x) * self.dt + self.diffusion() * np.random.normal(size=self.dim)



# Number of agents
num_agents = 9

# Number of dimensions (2 dimensions per agent)
dim = 2 * num_agents

# Initialize the state
# x = np.random.uniform(low=0.5, high=4000, size=(20,))
x = np.array([850,1420,900,1420,3870,580,3410,1710,3810,1350,1720,1300,2950,680,1600,2020,2690,1870],dtype='float64') # start near means
# x = np.random.rand(dim)

num_steps = 50000
time_step = 0.01

# Initialize an array to store the steps
steps = np.zeros((num_steps, dim))

# Create the GMM
gmm = GMM(means, covariances, weights) #Import means, covariances and weights from data

# Create the SDE
sde = SDE(2, time_step, gmm)  # Each agent is 2D

# Run the SDE for 1000 steps
for i in range(num_steps):
    # For each step, update each agent independently
    for j in range(num_agents):

        # Update the state of the current agent
        x[2*j:2*(j+1)] = sde.step(x[2*j:2*(j+1)])

    steps[i] = x

steps = steps.T


# Generate a grid of points
x1 = np.linspace(np.min(X_data[:,0]),np.max(X_data[:,0]),100)
y1 = np.linspace(np.min(X_data[:,1]),np.max(X_data[:,1]),100)
X, Y = np.meshgrid(x1, y1)

del x1, y1

# Calculate the GMM values at the grid points
Z = Gaussian_Mixture(X, Y)




plt.figure()
plt.clf()
# Plot the GMM contours
plt.contour(X, Y, Z, levels=50, cmap='viridis')

for i in range(num_agents):
    plt.plot(steps[2*i,:], steps[2*i+1,:])
plt.xlabel('Easting')
plt.ylabel('Northing')
plt.show()




plt.figure()
plt.plot(steps[0,:], steps[1,:])




from datetime import datetime

# Get the current date and time
now = datetime.now()

# Format date and time
date_str = now.strftime("%d%m%y")
time_str = now.strftime("%H%M")

# Format time_step to show only values after the decimal point
formatted_time_step = f"{time_step:.2f}".split('.')[1]

# Create the filename
filename = f"{date_str}_{time_str}_{num_agents}_{num_steps}_{formatted_time_step}.csv"

#Save to file and plot in MATLAB (check folder directory first to where it is being saved)
# date_time_num agents_steps_0.time step
np.savetxt(filename, steps.T, delimiter=',')
# np.savetxt('310724_1101_9_50000_01.csv',steps.T,delimiter = ',')

#%%

"""
In this code, dim is the total number of dimensions, which is twice the number of agents because each agent is represented by 2 dimensions. 
The state x and the steps steps are now arrays of size dim. 
After running the SDE, you can access the steps for a specific agent by indexing into the array. 
For example, steps[:, 0:2] gives all the steps for the first agent, steps[:, 2:4] gives all the steps for the second agent, and so on.

This gives you a 1x20 dimensional gmm - this code doesn't works!


"""

# import numpy as np
# from scipy.stats import multivariate_normal
# import matplotlib.pyplot as plt

# # Define the Gaussian Mixture Model
# class GMM:
#     def __init__(self, means, covariances, weights):
#         self.means = means
#         self.covariances = covariances
#         self.weights = weights

#     def potential(self, x):
#         return -np.log(sum(w * multivariate_normal(mean=m, cov=c).pdf(x) for m, c, w in zip(self.means, self.covariances, self.weights)))

# # Define the SDE
# class SDE:
#     def __init__(self, dim, dt, gmm):
#         self.dim = dim
#         self.dt = dt
#         self.gmm = gmm
#         self.eps = np.sqrt(np.finfo(float).eps)

#     def drift(self, x):
#         grad = np.zeros(self.dim)
#         for i in range(self.dim):
#             dx = np.zeros(self.dim)
#             dx[i] = self.eps
#             grad[i] = (self.gmm.potential(x + dx) - self.gmm.potential(x - dx)) / (2 * self.eps)
#         return -grad

#     def diffusion(self):
#         return np.sqrt(2 * self.dt)

#     def step(self, x):
#         return x + self.drift(x) * self.dt + self.diffusion() * np.random.normal(size=self.dim)

# # Number of agents
# num_agents = 10

# # Number of dimensions (2 dimensions per agent)
# dim = 2 * num_agents

# # Number of Gaussians
# num_gaussians = 8

# # Define the parameters of the GMM

# def generate_random_covariance(dim):
#     A = np.random.rand(dim, dim)
#     return np.dot(A, A.transpose())  # A*A^T is symmetric positive-definite


# mean = np.random.rand(dim)  # 2D mean
# # covariance = np.eye(dim)  # 2D covariance
# covariance = generate_random_covariance(dim)  # dim-dimensional covariance

# # Repeat the mean and covariance for each Gaussian
# means = [mean for _ in range(num_gaussians)]
# covariances = [covariance for _ in range(num_gaussians)]
# weights = [1/num_gaussians for _ in range(num_gaussians)]  # The weights should sum to 1

# # # Define the parameters of the GMM
# # means = [np.random.rand(dim) for _ in range(8)]
# # covariances = [np.eye(dim) for _ in range(8)]
# # weights = [1/8, 1/8, 1/8, 1/8,1/8, 1/8, 1/8, 1/8]

# # Create the GMM
# gmm = GMM(means, covariances, weights)

# # Create the SDE
# sde = SDE(dim, 0.01, gmm)

# # Initialize the state
# x = np.random.rand(dim)

# # Initialize an array to store the steps
# steps = np.zeros((1000, dim))

# # Run the SDE for 1000 steps
# for i in range(1000):
#     x = sde.step(x)
#     steps[i] = x

# # Now 'steps' is a numpy array that contains all the steps
# print(steps)

# steps = steps.T

# plt.figure()
# plt.clf()
# for i in range(num_agents):
#     plt.plot(steps[2*i,:], steps[2*i+1,:])
# plt.xlabel('Easting')
# plt.ylabel('Northing')

