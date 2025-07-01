
"""
Within this code we want to estimate the noise from data.

1. Import data. Require the paired co-ordinates (preferably in UTM format) and the corresponding time point. The CSV file should be checked before 
importing to ensure the time column has the format:

YYYY/MM/DD hh:mm:ss

2. Next we find the differenences in the time. This is followed by finding the differences in the locations. If there are n observations, then there are now n-1
observations. These two arrays are then combined.

3. Go through each column of the combined array of time and location differences and remove any columns where the time difference is above 120 minutes 
as this means there is a gap of 2 hours between observations. This can be tailored to your data depending on how often measurements are taken. Also remove any 
columns less than or equal to 0. This occurs if the column for time is not in the right format above or if estimating the noise for multiple animals and the date
and time is earlier than the prior animals reading.

4. Next, we split the array back up to have an array for the differences in location and another for the time. 

5. The next step is to estimate the matrix a = sigma * sigma.T 

6. Finally, apply Cholesky Decomposition to find an estimate for the matrix sigma in the form of a lower triangular matrix. Save the spydata for later use.

There is code at the bottom is the noise is to be estimated from generated data.

"""

import numpy as np
import pandas as pd


#%% Step 1: Import the Data

#Import Data - full data set
data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data From Dez/Data/Woodchester_UTM.csv')

x = data['Easting']
x = pd.Series.to_numpy(x)

y = data['Northing']
y = pd.Series.to_numpy(y, dtype='float64')

#Location data:
X = np.vstack((x,y)) 

data["DateTime"] = pd.to_datetime(data["DateTime"],)
Time = data['DateTime']

#Time data:
Time = pd.Series.to_numpy(Time)

del x,y #clean variable explorer

#%% Step 2: find the differences in time

delta_t = np.diff(Time).astype('timedelta64[m]') # difference in time
delta_t = (delta_t / np.timedelta64(1,'m')).astype(int) #Convert array of time from datetime to float (time difference calculated in m=minutes)
print(delta_t)


#%% Step 3: find the differences in location 
#X = distances row 1 = longitude/Easting, row 2 = latitude/Northing

delta_d = np.diff(X) # difference in distances

#%% Step 4: Combine and clean.

delta = np.vstack((delta_d,delta_t.T))

delta = np.delete(delta,(delta[2,:]>720), axis = 1) #Go through each column, if greater than 720 (i.e. a gap of more than 12 hours in the data) delete

delta = np.delete(delta,(delta[2,:]<=0), axis = 1) #Go through each column, reject each column that is less than or equal to 0

n = len(delta[0]) #total of observations left

[delta_d, delta_t] = np.array_split(delta,2, axis = 0) #Split the Data up to differences in location and differences in time

# del delta

#%% Steps 5 - 6: Estimate matrix a, then apply Cholesky decomposition
#% Change in Tau

A2 = np.zeros((2,2))

for i in range(n):
    A2 = A2 + (1/delta_t[i]) * np.outer(delta_d[:,i], delta_d[:,i]) 
A2 = 1/(n) * A2 #Applying Kramer's Moyal Formula to find an estimate for A

L2 = np.linalg.cholesky(A2)


#%% Reshape the data to get one value for sigma

delta_t1 = np.hstack((delta_t,delta_t))
delta_d1 = np.hstack((delta_d[0,:],delta_d[1,:]))

A1 = np.zeros((1))
for i in range(2*n):
    A1 = A1 + (1/delta_t1[:,i]) * np.outer(delta_d1[i], delta_d1[i]) 
A1 = 1/(2*n) * A1 #Applying Kramer's Moyal Formula to find an estimate for A

L1 = np.linalg.cholesky(A1)



#%% Using Generated Trajectory to Estimate Noise

"""
If using generated data, there is a slightly different format since there is no 'GMT Time'. 
Here, time is the step size, i.e. if step size is 0.1, the first entry is 0.1, the next is 0.2, and continues for all iterations.
Read in the data as normal, where the headings and stepsize might need to be added to the CSV file first.
Then, this code covers steps 1-4. The code does not need to be 'cleaned' as all entries are time entries one after the other. 
Run step 5-6 as usual.
"""

#read in data

data = pd.read_csv('/Users/jf01028/Library/CloudStorage/OneDrive-UniversityofSurrey/0.Data/Badger Data From Dez/EDMD/DATA/Badger_interp_175.csv')

x = data['Easting']
x = pd.Series.to_numpy(x)

y = data['Northing']
y = pd.Series.to_numpy(y, dtype='float64')

# Uncomment these lines if you have a time column in csv file
# Time = data['Time']
# Time = pd.Series.to_numpy(Time, dtype='float64')

# If no 'time column', adjust according to the step size
Time = np.zeros(len(y[0]))
for i in range(len(y[0])):
    Time[i] = (i+1) * 0.001

X = np.vstack((x,y))

del x,y, data 

#%% Calculate differences from original data and estimate diffusion

delta_t = np.diff(Time) #calculate differences in time

delta_d = np.diff(X) # difference in distances

n = len(delta_d[0])


delta_t1 = np.hstack((delta_t,delta_t))
delta_d1 = np.hstack((delta_d[0,:],delta_d[1,:]))

A1 = np.zeros((1))
for i in range(2*n):
    A1 = A1 + (1/delta_t1[:,i]) * np.outer(delta_d1[i], delta_d1[i]) 
A1 = 1/(2*n) * A1 #Applying Kramer's Moyal Formula to find an estimate for A

L1 = np.linalg.cholesky(A1)

#%% Estimate diffusion in 2x2 matrix

A2 = np.zeros((2,2))

for i in range(n):
    A2 = A2 + (1/delta_t[i]) * np.outer(delta_d[:,i], delta_d[:,i]) 
A2 = 1/(n) * A2 #Applying Kramer's Moyal Formula to find an estimate for A

L2 = np.linalg.cholesky(A2)

#%% To downsample data:
    
import random

data = np.vstack((X,Time)) #stack the coordinate data and the time data

n = len(data[0]) #find the length of the matrix

sel_cols = random.sample(range(n),500) #random selection of columns to keep

sel_cols.sort(reverse=False) #order lowest to highest

data_new = data[:, sel_cols] #new data set with reduced columns

del data, n, sel_cols

[distance_new, time_new] = np.array_split(data_new,2, axis = 0) 

delta_t_new1 = np.diff(time_new) #calculate differences in time

delta_d_new1 = np.diff(distance_new) # difference in distances

n = len(delta_d_new1[0])


