#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 13:31:25 2024

@author: jf01028
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from scipy.io import savemat
import pandas as pd

# Multivariate Gaussian PDF
def multivariate_gaussian_pdf(x, mean, cov):
    return multivariate_normal.pdf(x, mean, cov)

# Energy function E(x) = -log(p(x))
def energy_function(x, weights, means, covariances):
    p_x = sum(weights[k] * multivariate_gaussian_pdf(x, means[k], covariances[k]) for k in range(len(weights)))
    return -np.log(p_x)


# Responsibilities (posterior probabilities) w_k(x)
def compute_responsibilities(x, weights, means, covariances):
    K = len(weights)
    responsibilities = np.zeros(K)
    total = sum(weights[k] * multivariate_gaussian_pdf(x, means[k], covariances[k]) for k in range(K))
    total = max(total, 1e-10)  # Avoid zero division
    for k in range(K):
        responsibilities[k] = weights[k] * multivariate_gaussian_pdf(x, means[k], covariances[k]) / total
    return responsibilities

def first_derivative(x, weights, means, covariances):
    K = len(weights)
    N, D = x.shape
    derivatives = np.zeros((N, D))
    responsibilities = compute_responsibilities(x, weights, means, covariances)
    for i in range(N):
        for k in range(K):
            derivatives[i] += responsibilities[k] * np.dot(np.linalg.inv(covariances[k]), (x[i] - means[k]))
    return derivatives

# Gradient of the energy function
# def gradient_E(x, weights, means, covariances):
#     K = len(weights)
#     D = len(x)
#     grad = np.zeros(D)
#     responsibilities = compute_responsibilities(x, weights, means, covariances)
#     for k in range(K):
#         diff = x - means[k]
#         inv_cov_k = np.linalg.inv(covariances[k])
#         grad += responsibilities[k] * np.dot(inv_cov_k, diff)
#     return -grad  # Gradient should be negative for minimization

def gradient_E(x, weights, means, covariances):
    K = len(weights)
    D = len(x)
    grad = np.zeros(D)
    total_p_x = sum(weights[k] * multivariate_gaussian_pdf(x, means[k], covariances[k]) for k in range(K))
    total_p_x = max(total_p_x, 1e-10)  # Prevent division by zero

    for k in range(K):
        diff = x - means[k]
        inv_cov_k = np.linalg.inv(covariances[k])
        term = (weights[k] * multivariate_gaussian_pdf(x, means[k], covariances[k]) / total_p_x)
        grad += term * np.dot(inv_cov_k, diff)

    return -grad  # Negative sign for minimization


# Hessian of the energy function
def hessian_matrix(x, weights, means, covariances):
    K = len(weights)
    D = len(x)
    hessian = np.zeros((D, D))
    responsibilities = compute_responsibilities(x, weights, means, covariances)
    for k in range(K):
        diff = x - means[k]
        inv_cov_k = np.linalg.inv(covariances[k])
        term1 = inv_cov_k
        term2 = np.outer(np.dot(inv_cov_k, diff), np.dot(inv_cov_k, diff))
        hessian += responsibilities[k] * (term1 - term2)
    return hessian

# Use fsolve to find critical points
def find_critical_points(weights, means, covariances, initial_guesses):
    critical_points = []
    for guess in initial_guesses:
        critical_point, _, ier, _ = fsolve(gradient_E, guess, args=(weights, means, covariances), full_output=True,xtol=1e-10)
        if ier == 1:  # fsolve converged to a solution
            critical_points.append(critical_point)
    return np.array(critical_points)

# Compute Hessians at the critical points
def compute_hessians_at_critical_points(critical_points, weights, means, covariances):
    hessians = []
    for point in critical_points:
        hess = hessian_matrix(point, weights, means, covariances)
        hessians.append(hess)
    return np.array(hessians)

#%%

# Generate a grid of points for plot
# x_min, x_max = -1, 3
# y_min, y_max = -1, 3
x_min = np.min(X_data.T[0,:])
y_min = np.min(X_data.T[1,:])
x_max = np.max(X_data.T[0,:])
y_max = np.max(X_data.T[1,:])

x, y = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
z = np.zeros_like(x)


# Example GMM parameters
# weights = np.array([0.2, 0.3, 0.5])
# means = np.array([[0, 0], [1, 1], [2, 2]])
# covariances = np.array([np.eye(2), np.eye(2), np.eye(2)])

# Initial guesses for saddle points
# initial_guesses = [np.array([1500, 1500]), np.array([2000, 2000]), np.array([3500, 500])]
# initial_guesses = np.array([[1000,1500],[2000,1500],[1750,1500],[3500,500]])


x1 = np.linspace(0, round(x_max,-2), 20)
y1 = np.linspace(0,round(y_max,-2),20)
X_I, Y_I = np.meshgrid(x1,y1)
X_Flat = X_I.ravel()
Y_Flat = Y_I.ravel()
initial_guesses1 = np.vstack((X_Flat,Y_Flat)).T

del x1, y1, X_I, Y_I, X_Flat, Y_Flat

critical_points = find_critical_points(weights, means, covariances, initial_guesses1)
print("Critical Points:", critical_points)


# Step 1: Round the array to 2 decimal places
rounded_array_crit = np.round(critical_points, 2)

# Step 2: Sort the array (by rows)
# sorted_array_saddle_points = np.array([np.sort(row) for row in rounded_array_saddle_points])

# Step 3: Find unique rows (pairs)
df = pd.DataFrame(rounded_array_crit, columns=['col1', 'col2'])
unique_pairs_crit = df.drop_duplicates().values

# Compute Hessians at the critical points
hessians = compute_hessians_at_critical_points(unique_pairs_crit, weights, means, covariances)
print("Hessians at Critical Points:", hessians)

# Compute eigenvalues of the Hessians
eigenvalues = [np.linalg.eigvals(hess) for hess in hessians]
print("Eigenvalues of Hessians at Critical Points:", eigenvalues)

# Determine colors based on the signs of eigenvalues
colors = []
for eigvals in eigenvalues:
    if np.sign(eigvals[0]) != np.sign(eigvals[1]):
        colors.append('red')  # Saddle point (mountain pass)
    else:
        colors.append('black')  # Well

# Generate a grid of points
# x_min, x_max = -1, 3
# y_min, y_max = -1, 3
# x, y = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
# z = np.zeros_like(x)

# Compute the energy at each grid point
for i in range(x.shape[0]):
    for j in range(x.shape[1]):
        z[i, j] = energy_function(np.array([x[i, j], y[i, j]]), weights, means, covariances)

# Plot the contour of the energy surface
plt.contourf(x, y, z, levels=50, cmap='viridis')
plt.colorbar(label='Energy')

# Overlay the critical points with appropriate colors
for i, point in enumerate(unique_pairs_crit):
    plt.scatter(point[0], point[1], color=colors[i])

# Labels and title
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Energy Surface with Critical Points')
plt.legend()
plt.show()

roots = np.hstack((unique_pairs_crit,eigenvalues))

savemat('roots_means.mat', {'roots': roots})


# Save the dictionary to a .mat file
# savemat('hessians.mat', {f'matrix_{i+1}': hessians[i] for i in range(len(hessians))})

# savemat('saddle_points.mat', {'saddle_points': saddle_points})

# savemat('eigenvalues.mat', {'eigenvalues': eigenvalues})
