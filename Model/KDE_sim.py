
# nohup python3 <KDE_sim.py> output.txt &

#%% Import Spydata
import os
import sys
from spyder_kernels.utils.iofuncs import load_dictionary
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from datetime import datetime

# Get the current date and time
now = datetime.now()

# Format date and time
date_str = now.strftime("%d%m%y")
time_str = now.strftime("%H%M%S")


# Add the path to the custom module if necessary
sys.path.append('/user/HS502/jf01028/Documents/Simulations/d3s')

# Verify that the d3s module is available
try:
    import d3s
except ModuleNotFoundError:
    print("d3s module not found. Please ensure it is installed and accessible.")
    sys.exit(1)

# Path to your .spydata file
spydata_file_path = '/user/HS502/jf01028/Documents/Simulations/Data_M_C_W_KDE_withS.spydata'

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


#% Population class with the new potential

class Population(object):
    def __init__(self, KDE, Sigma):
        self.KDE = KDE #Generated Potential
        self.Sigma = Sigma #Estimated noise (Cholesky Decomposition = L)

    def b(self, x):     
        return -self.KDE.gradV(x) #negative gradient of the energy potential
                
    def sigma(self):
        return self.Sigma

# Define the SDE
class SDE:
    def __init__(self, dim, dt, KDE):
        self.dim = dim
        self.dt = dt
        self.KDE = KDE
        self.eps = np.sqrt(np.finfo(float).eps)
        self.bounds = [(0,4500),(0,2750)]

    def drift(self, x):
        grad = self.KDE.gradV(x)
        # print(f"gradV: {grad}")
        return -grad #negative gradient of the energy potential

    def diffusion(self):
        return 22 * np.sqrt( self.dt)
    
    def apply_boundary_conditions(self, x):
        for i in range(self.dim):
            min_bound, max_bound = self.bounds[i]
            if x[i] < min_bound:
                x[i] = 2 * min_bound - x[i]  # Reflect back into the range
            elif x[i] > max_bound:
                x[i] = 2 * max_bound - x[i]  # Reflect back into the range
        return x

    def step(self, x):
        new_x = x.T + self.drift(x[:,np.newaxis])[0,:] * self.dt + self.diffusion() * np.random.normal(size=self.dim)
        new_x = self.apply_boundary_conditions(new_x)
        return new_x

    # def step(self, x):
    #     return x.T + self.drift(x[:,np.newaxis])[0,:] * self.dt + self.diffusion() * np.random.normal(size=self.dim)




#% Initial conditions
# Number of agents
num_agents = 10

# # h = 0.05 # step size
# h = 0.1 # step size

# Number of dimensions (2 dimensions per agent)
dim = 2 * num_agents

num_steps = 1000000 # number of iterations

# Initialize the state
# x_kde1 = np.array([500, 1000,1000, 2000, 2000, 500, 3000, 2000,500, 1000,1000, 2000, 2000, 500, 3000, 2000, 600,1500,600,2000],dtype='float64') # 
x_kde = np.array([850,1420,900,1420,3870,580,3410,1710,3810,1350,1720,1300,2950,680,1600,2020,2690,1870,1915,915],dtype='float64') # start near means

# np.array([ 600,1500,600,1525],dtype='float64')
# np.random.uniform(low=0.5, high=4000, size=(dim,))
# x_kde = x_kde.reshape((1,2))

# # Initialize an array to store the steps
# steps = np.zeros((10000, dim))

time_step = 0.01

steps_kde = np.zeros((num_steps, dim))

# %
import time

start = time.time()

# Create the SDE
sde = SDE(2, time_step, KDE)  # type: ignore # Each agent is 2D

print('set-up completed')
sys.stdout.flush()

# Run the SDE for 1000 steps
for i in range(num_steps):
    # For each step, update each agent independently
    for j in range(num_agents):

        # Update the state of the current agent
        x_kde[2*j:2*(j+1)] = sde.step(x_kde[2*j:2*(j+1)])

    steps_kde[i] = x_kde

        
end = time.time()
print(end - start)

# import d3s.domain as domain

# xmin = np.min(X_data[:,0])
# ymin = np.min(X_data[:,1])
# xmax = np.max(X_data[:,0])
# ymax = np.max(X_data[:,1])
# bounds = np.array([[xmin, xmax],[ymin,ymax]])

# x1 = np.linspace(xmin,xmax,100)
# y1 = np.linspace(ymin,ymax,100)
# x1, y1 = np.meshgrid(x1, y1)


# #Set how many boxes you the potential to be visualised over. 
# #The larger the number the longer it takes, the smaller the number, the less accurate the potential.
# boxes = np.array([100, 100])

# Omega = domain.discretization(bounds, boxes) 

# PsiC = KDE.V(Omega.midpointGrid()) #potential
# PsiC_KDE1 = np.reshape(PsiC, (100, 100))

# Format time_step to show only values after the decimal point
formatted_time_step = f"{time_step:.2f}".split('.')[1]

# Specify the folder where you want to save the CSV file
folder_path = '/user/HS502/jf01028/Documents/Simulations/KDE'

# Ensure the folder exists
os.makedirs(folder_path, exist_ok=True)

# Create the filename
file_name = f"{date_str}_{time_str}_{num_agents}_{num_steps}_{formatted_time_step}.csv"

# Full path to the CSV file
file_path = os.path.join(folder_path, file_name)

#Save to file and plot in MATLAB (check folder directory first to where it is being saved)
# date_time_num agents_steps_0.time_step
np.savetxt(file_path, steps_kde, delimiter=',')
