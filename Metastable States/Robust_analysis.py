#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:19:21 2024

@author: jf01028
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import d3s.domain as domain
import d3s.observables as observables
import d3s.algorithms as algorithms

import time, sys

#%% Import Data
dataX = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/0.Woodchester/EDMD/Results/Woodchester_X.csv', header=None) #paste path file here
dataY = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/0.Woodchester/EDMD/Results/Woodchester_Y.csv',header=None) #paste path file here

site = 'WC_Total'

x = dataX[0]
x = pd.Series.to_numpy(x)

y = dataX[1]
y = pd.Series.to_numpy(y, dtype='float64')

X = np.vstack((x,y)) 

x = dataY[0]
x = pd.Series.to_numpy(x)

y = dataY[1]
y = pd.Series.to_numpy(y, dtype='float64')

Y = np.vstack((x,y)) 

del x,y, dataY, dataX #clean variable explorer

plt.scatter(*X)


ten_percent_batch_size = int(X.shape[1] * 0.1)

#%%

def read_csv_and_extract_every_10000(file_path):
    # Read the CSV file
    data = pd.read_csv(file_path)

    # Convert the first column to a numpy array
    first_column = data.iloc[:, 0].values

    # Extract the value after every 10,000 entries
    extracted_values = first_column[ten_percent_batch_size-1::ten_percent_batch_size]  # Start at index 9999, then step by 10,000

    return extracted_values

number_of_badgers = read_csv_and_extract_every_10000('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/0.Woodchester/EDMD/Results/Woodchester_interpolate.csv')
# print(extracted_values)

 #%%
print('import done')
sys.stdout.flush()
#% Set up EDMD
bounds = np.array([[round(np.min(X[0,:])-150,-1),round(np.max(X[0,:])+150,-1)],[round(np.min(X[1,:])-150,-1), round(np.max(X[1,:])+150,-1)]])

#%
box_x = np.round((bounds[0,1]-bounds[0,0])/100).astype(int)
box_y = np.round((bounds[1,1]-bounds[1,0])/100).astype(int)
boxes = np.array([box_x, box_y])

Omega = domain.discretization(bounds, boxes) 

C = Omega.midpointGrid()

# np.savetxt('C.csv',C, delimiter=',')

# choose observables
# psi = observables.monomials(10)
psi = observables.indicators(Omega)

evs = 40 # number of eigenvalues/eigenfunctions to be computed


print('setup done')
sys.stdout.flush()

PsiC = psi(Omega.midpointGrid()) # observables evaluated at midpoints of the grid
print('starting EDMD')
sys.stdout.flush()
#% Run EDMD
#% Streaming Version
import time
import scipy as sp

start = time.time()

def streaming_edmd(X, Y, psi, evs, operator='K', batch_size=1000):
    '''
    Streaming version of the conventional EDMD for the Koopman or Perron-Frobenius operator. 
    The matrices X and Y contain the input data.

    :param psi:       set of basis functions, see d3s.observables
    :param evs:       number of eigenvalues/eigenvectors
    :param operator:  'K' for Koopman or 'P' for Perron-Frobenius
    :param batch_size: size of the batches to process at a time
    :return:          eigenvalues d and corresponding eigenvectors V containing the coefficients for the eigenfunctions
    '''
    
    n_samples = X.shape[1]
    C_0 = np.zeros((psi(X[:, :batch_size]).shape[0], psi(X[:, :batch_size]).shape[0]))
    C_1 = np.zeros_like(C_0)
    
    # List to store A matrices from each batch
    # A_list = []
    d_list = []

    # Process data in batches
    for i in range(0, n_samples, batch_size):
        X_batch = X[:, i:i + batch_size]
        Y_batch = Y[:, i:i + batch_size]

        PsiX_batch = psi(X_batch)
        PsiY_batch = psi(Y_batch)

        C_0 += PsiX_batch @ PsiX_batch.T
        C_1 += PsiX_batch @ PsiY_batch.T
        
        # Compute and store A for this batch
        A_batch = sp.linalg.pinv(C_0) @ C_1
        # A_list.append(A_batch)
        d_batch , V = algorithms.sortEig(A_batch, evs)
        d_list.append(d_batch)

    if operator == 'P': 
        C_1 = C_1.T

    A = sp.linalg.pinv(C_0) @ C_1
    d, V = algorithms.sortEig(A, evs)
    
    return A, d, V, d_list


A, d, V, d_list = streaming_edmd(X, Y, psi, evs, operator='K', batch_size=ten_percent_batch_size)

end = time.time()
print(end - start)

np.savez(f'{site}_EDMD_Variables.npz',matrix1 = A, matrix2 = d , matrix3=V,matrix4 = d_list,matrix5 = number_of_badgers)

print('successfully saved variables A, d and V and number of badgers')
sys.stdout.flush()

#%

#%
# import numpy as np

# Load the .npz file
# data = np.load('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data From Dez/EDMD/Robust/WC_2019_EDMD_Variables.npz')

# # Access the matrices and other variables
# A = data['matrix1']
# d = data['matrix2']
# V = data['matrix3']
# d_list = data['matrix4']
# number_of_badgers = data['matrix5']


#%

import numpy as np
import matplotlib.pyplot as plt
import os

def plot_eigenvalues(d_list, output_folder):
    """
    Plot eigenvalues for each batch.

    :param eigenvalues_list: List of eigenvalues computed at each batch.
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    num_batches = len(d_list)
    
    n = 21 # Display n-1 eigenvalues
    
    for i in range(num_batches):
        eigenvalues = d_list[i].real
        
        # Calculate the differences between consecutive eigenvalues
        EV_diff = np.zeros((len(eigenvalues) - 1, 1))
        for j in range(len(eigenvalues) - 1):
            EV_diff[j] = abs(eigenvalues[j] - eigenvalues[j + 1])
        
        # Create subplots with shared x-axis
        fig, ax1 = plt.subplots(figsize=(8, 4))
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis
        
        # Plot real eigenvalues on the first axis
        ax1.plot(eigenvalues[0:n], 'g-o', label='Real Part of Eigenvalues')
        
        # Plot the differences on the second axis
        ax2.plot(EV_diff[0:n], 'b-o', label='Eigenvalue Differences')
        
        # Set axis labels and grid
        ax1.set_xlabel('Eigenvalue Index')
        ax1.set_ylabel('Real Part of Eigenvalues', color='g')
        ax2.set_ylabel('Eigenvalue Differences', color='b')
        
        # Set grid and ticks
        ax1.grid(which='major', axis='both', linestyle='--')
        ax1.set_xticks(range(n))  # Set x-ticks as integers
        ax1.set_xticklabels(range(n))  # Set x-tick labels as integers
        
        # Show plot with titles and legend
        plt.title(f'{site}_Eigenvalues and Differences - Batch {i+1}/{num_batches}')
        fig.tight_layout()  # Adjust layout to prevent overlap
        # plt.show()

        # Save the plot to the specified folder
        plot_filename = os.path.join(output_folder, f'{site}_eigenvalues_batch_{i+1}.png')
        plt.savefig(plot_filename)
        
        # Close the plot to free memory
        plt.close(fig)

# Assuming d_list is already defined and contains the eigenvalues
output_folder = '/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/0.Woodchester/EDMD/Robust/Total'
plot_eigenvalues(d_list,output_folder)





#%%

# WC_SPEC = np.array([ 2,  2,  2,  3,  6,  7,  7,  9,  9,  7,  7, 10, 10,  7,  8,  8,  8,  8,  8,  8,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
#        10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13,  8,  8,  8,  8, 8,  8,  8,  8,  8,  8,  8,  8,  8])



# # Create subplots with shared x-axis
# fig, ax1 = plt.subplots(figsize=(8, 4))
# ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis

# # Plot real eigenvalues on the first axis
# ax1.plot(WC_SPEC, 'g-o', label='Real Part of Eigenvalues')

# # Plot the differences on the second axis
# ax2.plot(number_of_badgers, 'b-o', label='Number of Badgers')

# # Set axis labels and grid
# ax1.set_xlabel('Batch Number')
# ax1.set_ylabel('Spectral Gap - metastable states', color='g')
# ax2.set_ylabel('Count of Badgers', color='b')

# # Set grid and ticks
# ax1.grid(which='major', axis='both', linestyle='--')
# # plt.xticks(range(len(eigenvalues)))  # Adjust ticks as needed

# # Show plot with titles and legend
# plt.title('Metastable States and Badgers')
# fig.tight_layout()  # Adjust layout to prevent overlap
# plt.show()



