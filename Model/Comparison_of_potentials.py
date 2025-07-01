#%% Import Spydata and define classes for SDE
import os
import sys
# from spyder_kernels.utils.iofuncs import load_dictionary
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import multivariate_normal
from datetime import datetime
# from numba import jit
import pickle
from scipy.optimize import minimize


# Load the dictionary from the file
with open('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/2.CODING/dynamical-systems/Comparison of Potentials/variables.pkl', 'rb') as f:
    variables = pickle.load(f)

# Access the variables from the dictionary
covariances = variables["covariances"]
gmm = variables["gmm"]
KDE = variables["KDE"]
labels = variables["labels"]
means = variables["means"]
num_agents = variables["num_agents"]
num_steps = variables["num_steps"]
S = variables["S"]
test = variables["test"]
weights = variables["weights"]
X = variables["X"]
X_data = variables["X_data"]


def gmm_energy_potential(x, means, covariances, weights):
    return -(22^2)/2 * np.log(sum(w * multivariate_normal(mean=m, cov=c).pdf(x) for m, c, w in zip(means, covariances, weights)))

def kde_energy_potential(x, KDE):
    return KDE.V(x.reshape(2, 1))[:,0]

# Test the potentials
x = np.array([0.5, 0.5])

print('GMM Potential:', gmm_energy_potential(x,means,covariances,weights))
print('KDE Potential:', kde_energy_potential(x,KDE)[0])

# Define the objective function
def objective_function(c,x,means,covariances,weights,KDE):
    gmm_pot = gmm_energy_potential(x, means, covariances, weights)
    kde_pot = kde_energy_potential(x, KDE) + c
    return abs(gmm_pot - kde_pot)

# Fixed point x
x = np.array([1000.0, 1000.0])

# Initial guess for c
c0 = 0.0

# Define bounds for c (if any)
bounds = [(None, None)]  # No bounds in this example

# Perform the optimization
result = minimize(objective_function, c0, args=(x, means, covariances, weights, KDE), bounds=bounds, method='SLSQP')

# Display results
optimal_c = result.x[0]
print('Optimal solution for c:', optimal_c)
print('Objective function value at optimal solution:', result.fun)

# %%

# Define the center point
center_point = means[0]

# Define the range around the center point
range_x = 50  # Range in the x direction
range_y = 50  # Range in the y direction

# Define the number of points in each direction
num_points = 100

# Create the domain around the center point
x_range = np.linspace(center_point[0] - range_x / 2, center_point[0] + range_x / 2, num_points)
y_range = np.linspace(center_point[1] - range_y / 2, center_point[1] + range_y / 2, num_points)
X, Y = np.meshgrid(x_range, y_range)
grid_points = np.c_[X.ravel(), Y.ravel()]

# Plot the domain and the center point
plt.figure(figsize=(8, 8))
plt.scatter(grid_points[:, 0], grid_points[:, 1], s=1, label='Grid Points')
plt.scatter(center_point[0], center_point[1], color='red', label='Center Point')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.title('Box Domain with Center Point')
plt.show()


#%%



# Define the domain
x_range = np.linspace(0, 50, 100)
y_range = np.linspace(0, 50, 100)
X, Y = np.meshgrid(x_range, y_range)
grid_points = np.c_[X.ravel(), Y.ravel()]

# Evaluate the potentials on the grid
gmm_potentials = np.array([gmm_energy_potential(point, means, covariances, weights) for point in grid_points])
kde_potentials = np.array([kde_energy_potential(point, KDE) for point in grid_points])

# Compute the difference
differences = gmm_potentials - kde_potentials

# Integrate over the domain (sum the differences weighted by the area of each grid cell)
dx = x_range[1] - x_range[0]
dy = y_range[1] - y_range[0]
area = dx * dy
integral = np.sum(differences**2) * area

print('Integral of the squared difference over the domain:', integral)
# %%
