#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimate the Diffusion term for the population without splitting the file into individuals/years/months

Must have columns named 'Year' for the year, 'Month' for the month and 'Badger_Ind' for badger identity. These can be changed based on own file.

Must also have 'Easting' and 'Northing' column, and 'DateTime' column in the following format: YYYY/MM/DD hh:mm:ss

The code estimates a matrix for diffusion and computes the Cholesky Decomposition and also a single value for the diffusion, c (\sigma(x)=cI)

Outputs into one file.

"""

import numpy as np
import pandas as pd

data = pd.read_csv('file.csv') # import data


newlist = []

Year_list = data['Year'].unique().tolist()


for year in Year_list:
    
    # Filter for a specific year
    data_y = data[data['Year']==year]
    
    Month_list = data_y['Month'].unique().tolist()
    
    for month in Month_list:
        
        data_m = data_y[data_y['Month']==month]

        ID_list = data_m['Badger_Ind'].unique().tolist() #column named with the ID of the badgers

        for ID in ID_list:
            innerlist = []
        
            data1 = data_m[data_m['Badger_Ind']==ID]
            
            animal = data1['Badger_Ind'] #the column that is the same for each csv file
            
            innerlist.append(animal.iloc[0])
            
            # if doing multiple sites at once, uncomment the two lines below:
            
        #     # site = data['Site']
            
        #     # innerlist.append(site[0])
            
            Year = data1['Year']
            Month = data1['Month']
            Sex = data1['Sex']
            
            innerlist.append(Year.iloc[0])
            innerlist.append(Month.iloc[0])
            innerlist.append(Sex.iloc[0])
            
            newlist.append(innerlist)
            
            x = data1['Easting']
            x = pd.Series.to_numpy(x)

            y = data1['Northing']
            y = pd.Series.to_numpy(y, dtype='float64')

            #Location data:
            X = np.vstack((x,y)) 

            data1["DateTime"] = pd.to_datetime(data1["DateTime"],)
            Time = data1['DateTime']

            #Time data:
            Time = pd.Series.to_numpy(Time)

            del x,y #clean variable explorer
    
            innerlist.append(len(X[0]))
    
            #Step 2: find the differences in time
        
            delta_t = np.diff(Time).astype('timedelta64[m]') # difference in time
            delta_t = (delta_t / np.timedelta64(1,'m')).astype(int) #Convert array of time from datetime to float (time difference calculated in m=minutes)
            # print(delta_t)

            # Step 3: find the differences in location 
            #X = distances row 1 = longitude/Easting, row 2 = latitude/Northing

            delta_d = np.diff(X) # difference in distances

            #% Step 4: Combine and clean.
        
            delta = np.vstack((delta_d,delta_t.T))
        
            delta = np.delete(delta,(delta[2,:]>720), axis = 1) #Go through each column, if greater than 720 (i.e. a gap of more than 12 hours in the data) delete
        
            delta = np.delete(delta,(delta[2,:]<=0), axis = 1) #Go through each column, reject each column that is less than or equal to 0
        
            n = len(delta[0]) #total of observations left
        
            [delta_d, delta_t] = np.array_split(delta,2, axis = 0) #Split the Data up to differences in location and differences in time
        
            innerlist.append(n)

            # Steps 5: Estimate matrix a, then apply Cholesky decomposition
            
            A = np.zeros((2,2))
        
            for i in range(n):
                A = A + (1/delta_t[:,i]) * np.outer(delta_d[:,i], delta_d[:,i]) 
            A = 1/(n) * A #Applying Kramer's Moyal Formula to find an estimate for A
        
            if n > 1:
                L = np.linalg.cholesky(A) 
            else:
                L = np.zeros((2,2))
    
            
            innerlist.append(A[0,0])
            innerlist.append(A[0,1])
            innerlist.append(A[1,0])
            innerlist.append(A[1,1])
            
            F = np.linalg.norm(L,ord='fro')
                    
            innerlist.append(L[0,0])
            innerlist.append(L[0,1])
            innerlist.append(L[1,0])
            innerlist.append(L[1,1])
            
            innerlist.append(F)
            
            # Step 5b: Estimate single value for diffusion c (i.e. \sigma(x) = cI)
            
            delta_t = np.hstack((delta_t,delta_t))
            delta_d = np.hstack((delta_d[0,:],delta_d[1,:]))
        
            A = np.zeros((1))
            for i in range(2*n):
                A = A + (1/delta_t[:,i]) * np.outer(delta_d[i], delta_d[i]) 
            A = 1/(2*n) * A #Applying Kramer's Moyal Formula to find an estimate for A
        
            if n > 1:
                L = np.linalg.cholesky(A) # i.e. sqrt(A)
            else:
                L = np.zeros((2,2))
            
            innerlist.append(np.take(A ,0))
            innerlist.append(np.take(L,0))
            
            sigma = L*np.identity((2))
            
            F1 = np.linalg.norm(sigma,ord='fro') 
            
            innerlist.append(F1) 
            
            # if interested in season, this can be added in too:
            # season = data['Season']
            # innerlist.append(season[0])
    

column_names = ["Animal","Year","Month", "Sex","No Data Points", "Used points","A11","A12","A21","A22","L11","L12","L21","L22","Lfro_matrix","c_squared","c","Lfro_single"]

df = pd.DataFrame(newlist, columns=column_names)
df.to_csv('Site.csv', index=False, header=True)
