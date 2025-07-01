
# nohup python3 <KDE_sim.py> output.txt &

#%% Import Spydata
import os
import sys
from spyder_kernels.utils.iofuncs import load_dictionary
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import multivariate_normal
from datetime import datetime

# Get the current date and time
now = datetime.now()

# Format date and time
date_str = now.strftime("%d%m%y")
time_str = now.strftime("%H%M%S")


# Add the path to the custom module if necessary
sys.path.append('/Simulations/d3s')

# Verify that the d3s module is available
try:
    import d3s
except ModuleNotFoundError:
    print("d3s module not found. Please ensure it is installed and accessible.")
    sys.exit(1)

# Path to your .spydata file
spydata_file_path = '/Simulations/Data_M_C_W_KDE_withS.spydata'

# Load the .spydata file
def load_spydata(file_path):
    if os.path.exists(file_path):
        data, error = load_dictionary(file_path)
        if error:
            print(f"Error loading spydata file: {error}")
        else:
            return data
    else:
        print(f"File not found: {file_path}")
        return None

# Load the data
data = load_spydata(spydata_file_path)

# Print the data
if data:
    globals().update(data)
    print("Loaded variables:")
    for key in data.keys():
        print(key)
# %

print('spydata imported')
sys.stdout.flush()

#%
#% Population class with the new potential


#To plot the contours:
def Gaussian_Mixture(x, y):
    return np.sum([weights[i] * multivariate_normal.pdf(np.dstack([x, y]), mean=means[i], cov=covariances[i]) for i in range(len(weights))], axis=0)# type: ignore


class GMM:
    def __init__(self, means, covariances, weights):
        self.means = means
        self.covariances = covariances
        self.weights = weights

    def potential(self, x):
        return -np.log(sum(w * multivariate_normal(mean=m, cov=c).pdf(x) for m, c, w in zip(self.means, self.covariances, self.weights)))
    
    
# Define the SDE
class SDE: # type: ignore
    def __init__(self, dim, dt, gmm):
        self.dim = dim #dimension
        self.dt = dt #time step
        self.gmm = gmm
        self.eps = np.sqrt(np.finfo(float).eps)
        self.bounds = [(0,4500),(0,2750)]

    def multivariate_gaussian_pdf(self, xv):
        return multivariate_normal.pdf(x, means, covariances) # type: ignore

    # Responsibilities (posterior probabilities) w_k(x)
    def compute_responsibilities(self, x):
        K = len(weights) # type: ignore
        responsibilities = np.zeros(K)
        total = sum(weights[k] * self.multivariate_gaussian_pdf(x, means[k], covariances[k]) for k in range(K)) # type: ignore
        for k in range(K):
            responsibilities[k] = weights[k] * self.multivariate_gaussian_pdf(x, means[k], covariances[k]) / total # type: ignore
        return responsibilities

    def gradient_E(self, x): #works when only 1 agent
        K = len(weights) # type: ignore
        D = len(x)
        grad = np.zeros(D)
        responsibilities = self.compute_responsibilities(x, weights, means, covariances) # type: ignore
        for k in range(K):
            diff = x - means[k] # type: ignore
            inv_cov_k = np.linalg.inv(covariances[k]) # type: ignore
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
    
    def apply_boundary_conditions(self, x):
        for i in range(self.dim):
            min_bound, max_bound = self.bounds[i]
            if x[i] < min_bound:
                x[i] = 2 * min_bound - x[i]  # Reflect back into the range
            elif x[i] > max_bound:
                x[i] = 2 * max_bound - x[i]  # Reflect back into the range
        return x

    def step(self, x):
        new_x = x + self.drift(x) * self.dt + self.diffusion() * np.random.normal(size=self.dim)
        new_x = self.apply_boundary_conditions(new_x)
        return new_x

#%


#% Initial conditions
# Number of agents
num_agents = 10


# Number of dimensions (2 dimensions per agent)
dim = 2 * num_agents

num_steps = 100000 # number of iterations
time_step = 0.01

# Initialize the state
# x_kde1 = np.array([500, 1000,1000, 2000, 2000, 500, 3000, 2000,500, 1000,1000, 2000, 2000, 500, 3000, 2000, 600,1500,600,2000],dtype='float64') # 
x = np.array([850,1420,900,1420,3870,580,3410,1710,3810,1350,1720,1300,2950,680,1600,2020,2690,1870,1915,915],dtype='float64') # start near means

# x = np.array([4450,1500],dtype='float64')

# np.array([ 600,1500,600,1525],dtype='float64')
# np.random.uniform(low=0.5, high=4000, size=(dim,))
# x_kde = x_kde.reshape((1,2))

# # Initialize an array to store the steps
# steps = np.zeros((10000, dim))


# Initialize an array to store the steps
steps = np.zeros((num_steps, dim))

# Create the GMM
gmm = GMM(means, covariances, weights) # type: ignore #Import means, covariances and weights from data

# Create the SDE
sde = SDE(2, time_step, gmm)  # Each agent is 2D

print('set-up completed')
sys.stdout.flush()

# %
import time

start = time.time()

# Run the SDE for 1000 steps
for i in range(num_steps):
    # For each step, update each agent independently
    for j in range(num_agents):

        # Update the state of the current agent
        x[2*j:2*(j+1)] = sde.step(x[2*j:2*(j+1)])

    steps[i] = x

steps = steps.T

end = time.time()
print(end - start)

# #%%

# # Generate a grid of points
# x1 = np.linspace(np.min(X_data[:,0]),np.max(X_data[:,0]),100) # type: ignore
# y1 = np.linspace(np.min(X_data[:,1]),np.max(X_data[:,1]),100) # type: ignore
# X, Y = np.meshgrid(x1, y1)

# del x1, y1

# # Calculate the GMM values at the grid points
# Z = Gaussian_Mixture(X, Y)

# plt.figure()
# plt.clf()
# # Plot the GMM contours
# plt.contour(X, Y, Z, levels=50, cmap='viridis')

# for i in range(num_agents):
#     plt.plot(steps[2*i,:], steps[2*i+1,:])
# plt.xlabel('Easting')
# plt.ylabel('Northing')
# plt.show()


# plt.figure()
# plt.plot(steps[0,:], steps[1,:])

#%

# Format time_step to show only values after the decimal point
formatted_time_step = f"{time_step:.2f}".split('.')[1]

# Specify the folder where you want to save the CSV file
folder_path = '/user/HS502/jf01028/Documents/Simulations/GMM'

# Ensure the folder exists
os.makedirs(folder_path, exist_ok=True)

# Create the filename
file_name = f"GMM_{date_str}_{time_str}_{num_agents}_{num_steps}_{formatted_time_step}.csv"

# Full path to the CSV file
file_path = os.path.join(folder_path, file_name)

#Save to file and plot in MATLAB (check folder directory first to where it is being saved)
# date_time_num agents_steps_0.time_step
np.savetxt(file_path,steps.T, delimiter=',')
