#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:37:35 2024

@author: jf01028
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

import scipy
import scipy.spatial.distance

X = pd.read_csv('.csv',header=None)


X = X.to_numpy()

# Number of agents
n = X.shape[1] // 2

# Plot the coordinates
plt.figure(figsize=(10, 10))
for i in range(n):
    plt.scatter(X[:, 2*i], X[:, 2*i + 1], label=f'Agent {i+1}')

plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Agent Trajectories')
# plt.legend()
plt.show()

#%%

# Function to find the first time step where any agents are within 1 unit of each other
def find_first_close_encounter(X, n, threshold=5):
    for t in range(X.shape[0]):
        # Extract the coordinates at time step t
        coords = X[t, :].reshape(n, 2)
        
        # Compute pairwise distances
        dists = scipy.spatial.distance.pdist(coords)
        
        # Check if any distances are less than or equal to the threshold
        if np.any(dists <= threshold):
            # Find the indices of the agents involved
            dist_matrix = scipy.spatial.distance.squareform(dists)
            close_pairs = np.where(dist_matrix <= threshold)
            agent_pairs = [(i, j) for i, j in zip(close_pairs[0], close_pairs[1]) if i < j]
            return t, coords, agent_pairs

    return None, None, None

# Find the first close encounter
time_step, coords, agent_pairs = find_first_close_encounter(X, n)

if time_step is not None:
    print(f"The first close encounter occurs at time step {time_step}.")
    print("Coordinates of agents at this time step:")
    print(coords)
    print("Agent pairs within 1 unit of each other:")
    for pair in agent_pairs:
        print(f"Agent {pair[0]+1} and Agent {pair[1]+1}")
else:
    print("No close encounters found within the given threshold.")



#%% Transition of S to I

"""
The next step is to calculate the transition from S to I. 

A matrix is calculated to determine what category each badger is in. This is done for each time point. Then, if the badger is within a certain distance (r)
than the badger will make the transition from S to I.

S is susceptible and I is infectious.

"""

nIter = X.shape[0]  # Number of time points

Y = np.zeros((n, nIter), dtype = str)

Y[:,0] = np.full((1, n), 'S', dtype=str)
Y[0,0] = 'I' #First badger is infected and infectious

r = 1
# Iterate over each time step
for t in range(nIter - 1):
    Y[:, t + 1] = Y[:, t]
    
    # Extract the coordinates at time step t
    coords = X[t, :].reshape(n, 2)
    
    # Compute pairwise distances
    dist_matrix = scipy.spatial.distance.cdist(coords, coords)
    
    for i in range(n):
        if Y[i, t] == 'I':
            # Find indices of agents within distance r
            ind = np.where((dist_matrix[i, :] < r) & (dist_matrix[i, :] > 0))[0]
            for y in ind:
                p = np.random.rand()
                if p > 0.75:
                    Y[y, t + 1] = 'I'

# Initialize arrays to hold the counts of S and I for each time point
count_S = np.zeros(Y.shape[1])
count_I = np.zeros(Y.shape[1])

for t in range(Y.shape[1]):
    count_S[t] = np.sum(Y[:, t] == 'S')  # Count the number of 'S' in the column
    count_I[t] = np.sum(Y[:, t] == 'I')  # Count the number of 'I' in the column

time_points = np.arange(Y.shape[1])

# Plot the counts of S and I over time
plt.figure(figsize=(10, 6))
plt.plot(time_points, count_S, label='S (Susceptible)', color='b')
plt.plot(time_points, count_I, label='I (Infected)', color='r')

plt.xlabel('Time Points')
plt.ylabel('Counts')
plt.title('Counts of Susceptible and Infected Agents Over Time')
plt.legend()
# Set y-ticks to show every 2 units
plt.yticks(np.arange(0, max(count_S.max(), count_I.max()) + 2, 2))
plt.grid(True)
plt.show()



#%%


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance

# Define the folder containing the files
folder_path = 'First5'

# List all files in the folder
files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

# Initialize lists to hold the counts of S and I for each file
all_counts_S = []
all_counts_I = []

# Iterate over each file
for file in files:
    # Read the file into a DataFrame
    file_path = os.path.join(folder_path, file)
    X = pd.read_csv(file_path, header=None).to_numpy()

    # Number of agents
    n = X.shape[1] // 2

    # Initialize the infection status array
    nIter = X.shape[0]  # Number of time points
    Y = np.zeros((n, nIter), dtype=str)
    Y[:, 0] = np.full((1, n), 'S', dtype=str)
    Y[0, 0] = 'I'  # First agent is infected and infectious

    r = 1  # Infection radius

    # Iterate over each time step
    for t in range(nIter - 1):
        Y[:, t + 1] = Y[:, t]

        # Extract the coordinates at time step t
        coords = X[t, :].reshape(n, 2)

        # Compute pairwise distances
        dist_matrix = scipy.spatial.distance.cdist(coords, coords)

        for i in range(n):
            if Y[i, t] == 'I':
                # Find indices of agents within distance r
                ind = np.where((dist_matrix[i, :] < r) & (dist_matrix[i, :] > 0))[0]
                for y in ind:
                    p = np.random.rand()
                    if p > 0.9:
                        Y[y, t + 1] = 'I'

    # Count the number of S and I at each time point
    count_S = np.zeros(Y.shape[1])
    count_I = np.zeros(Y.shape[1])

    for t in range(Y.shape[1]):
        count_S[t] = np.sum(Y[:, t] == 'S')  # Count the number of 'S' in the column
        count_I[t] = np.sum(Y[:, t] == 'I')  # Count the number of 'I' in the column

    all_counts_S.append(count_S)
    all_counts_I.append(count_I)

# Find the maximum length of the counts
max_length = max(len(count) for count in all_counts_S)

# Pad the counts with zeros to make them the same length
all_counts_S = np.array([np.pad(count, (0, max_length - len(count)), 'constant') for count in all_counts_S])
all_counts_I = np.array([np.pad(count, (0, max_length - len(count)), 'constant') for count in all_counts_I])

# Plot the counts of S and I over time for each file
time_points = np.arange(max_length)

plt.figure(figsize=(10, 6))
for i in range(len(files)):
    plt.plot(time_points, all_counts_I[i], label=f'File {i+1}')

plt.xlabel('Time Points')
plt.ylabel('Number of Infected Agents')
plt.title('Number of Infected Agents Over Time')
plt.legend()
plt.grid(True)
plt.show()

#%%

# Select every 43800th data point (12 time points)
step = 43800
selected_columns = all_counts_I[:, ::step]  # Slice every 43800th column

# Calculate the number of infections for each selected time point
num_infections = np.sum(selected_columns, axis=0)

print(f"Number of infections at each 43800th time point: {num_infections}")



#%% Power Analysis

from statsmodels.stats.power import FTestAnovaPower

# Define parameters
effect_size = 0.4  # Large effect size
alpha = 0.05        # Significance level
power = 0.8         # Desired power
num_groups = 12     # Number of time points

# Perform power analysis
analysis = FTestAnovaPower()
sample_size = analysis.solve_power(effect_size=effect_size, alpha=alpha, power=power, k_groups=num_groups)

print(f"Required sample size per group: {sample_size:.2f}")
print(f"Total samples required: {sample_size * num_groups:.2f}")

#  For a one way anova - investigating if the average number of infected agents does not change significantly over time
# Required sample size per group: 115.19
# Total samples required: 1382.30

#%% Animation

"""
Then to view the transitions, we create an animation.
"""

fig = plt.figure(1)
fig.clf()
ims = []
for t in range(5000):
    # c = ['#0000ff' if Y[i, t] == 'S' else '#ff0000' for i in range(n)]
    im = plt.scatter(X[t,::2], X[t,1::2], s=20,color='red')
    
    im.axes.set_xlim(-2, 4500)
    im.axes.set_ylim(-1, 2750)
    ims.append([im])
    
# plt.axvline(x = 0, ymin = 0, ymax = 3, ls = '--', label = 'Border Line') #The boarder line (x=0)
# plt.axhline(y = 0, xmin = -4, xmax = 3, ls = '--') #The boarder line (y=0)

#Double Well
# plt.plot(1,0,'go', label = 'P.O.A., i.e. Sett')
# plt.plot(-1,0,'go')

#Quadruple Well
# plt.plot(1,-1,'go', label = 'P.O.A., i.e. Sett')
# plt.plot(-1,-1,'go')
# plt.plot(-1,1,'go')
# plt.plot(1,1,'go')

# plt.legend(loc='lower left')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

anim = animation.ArtistAnimation(fig, ims, interval = 10, blit = True)
anim.save('badgers_infection.mp4', writer='ffmpeg')
