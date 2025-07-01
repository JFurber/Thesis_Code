#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code investigates multiple badgers in a double-well energy potential with a simple susceptible-infectious model on top.

Used in the chapter  'Modelling Strategies', under the section ' Infection Model'.

The first half generates a trajectory as seen in 'Badgers.py' and 'EcologyModel.py' in an energy potential.

The second half add an SI model.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import scipy
import scipy.spatial


from datetime import datetime

date = datetime.now().strftime("%Y_%m_%d")


#%%
class Population(object):
    def __init__(self, A, Sigma, alpha, beta):
        self.A = A #behaviour matrix
        self.n = A.shape[0]
        self.Sigma = Sigma
        self.alpha = alpha #weight of attraction
        self.beta = beta #weight of repulsion

    def Attraction(self, x, i, j): 
        return -2*(x[2*i:2*i+2] - x[2*j:2*j+2])
    
    def Repulsion(self, x, i, j):
        return np.array([ 2*(x[2*i] - x[2*j])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2, 2*(x[2*i+1] - x[2*j+1])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2])
        

    def b(self, x):
        y = np.zeros_like(x) # Returns an array of given shape and type as the given array (i.e. x), with zeros
        for i in range(0, self.n):
            y[2*i:2*i+2] = np.array([-4*x[2*i]**3 + 4*x[2*i] , -2*x[2*i+1]]) #Double well potential
            # y[2*i:2*i+2] = np.array([-4*x[2*i]**3 + 4*x[2*i] , -4*x[2*i+1]**3 + 4*x[2*i+1]]) #Quadruple well potential
            for j in range(0, self.n):
                if self.A[i][j] == 1:
                    y[2*i:2*i+2] += self.alpha * self.Attraction(x, i, j)
                elif self.A[i][j] == -1:
                    y[2*i:2*i+2] += self.beta * self.Repulsion(x, i, j)
        return y
        
    def sigma(self):
        return self.Sigma

#%% number of badgers
n = 10 # How many badgers are in the system
f = 4
m = n - f

# %% Random generated A - behaviour matrix
A=2*np.random.randint(0,2,size =(n,n))-1
np.fill_diagonal(A, 0)
print(A)

#%% Generation of Noise - under the assumption the first f badgers are female and the rest are male
"""
This is not currently random, but can be changed to be random if need be, depending on the amount of badgers.
"""

f_noise = 0.5 #Noise Term for females
m_noise = 0.9 #Noise Term for males

F = np.full((2*f,2*f),f_noise)
M = np.full((2*m,2*m),m_noise)
a1 = np.diag(F)
a2 = np.diag(M)
a3 = np.concatenate((a1,a2))
Sigma = np.diag(a3)
del F, M, a1, a2, a3

# %% Population 
P = Population(A, Sigma ,alpha = 0.01, beta = 0.01)

# Initial conditions
h = 1e-3 # step size
d = 2*n # State space
nIter = 5000 # number of iterations

X = np.zeros((d, nIter))
X[:,0] = 4*np.random.rand(d)-2

#%% Euler-Maruyama
for t in range(1, nIter):
    X[:,t] = X[:,t-1] + P.b(X[:,t-1])*h + P.sigma()@np.random.randn(d)*np.sqrt(h)

# Plot of Results
plt.figure(1)
plt.clf()
for i in range(n):
    plt.plot(X[2*i,:], X[2*i+1,:])
    plt.plot(X[2*i,0], X[2*i+1,0],'x')


"""
From here we focus on generating an SI Model. 

"""
#%% Compute Distances

"""
The first step is calculating the Euclidean distance between each badger.

This is done at each time step and the distances are saved in a matrix, which are saved in a list. 
"""
    
D = []
for t in range(0,nIter):
    x1=np.zeros((n,2))
    for k in range(3):
            x1[:,0] = X[::2,t]
            x1[:,1] = X[1::2,t]
    Di = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(x1)) #distance is seen in a matrix
    D.append([Di])
    del Di

#%% Transition of S to I

"""
The next step is to calculate the transition from S to I. 

A matrix is calculated to determine what category each badger is in. This is done for each time point. Then, if the badger is within a certain distance (r)
than the badger will make the transition from S to I.

S is susceptible and I is infectious.

"""

Y = np.zeros((n, nIter), dtype = str)

Y[:,0] = np.full((1, n), 'S', dtype=str)
Y[0,0] = 'I' #First badger is infected and infectious

# Y[: , 0] = np.array(['S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'I'])

r = 0.25 #Distance to catch

for t in range (0,nIter-1):
    Y[:,t+1] = Y[:,t]
    
    for i in range (0,n):
        if Y[i,t] == 'I':
            ind = np.where(D[t][0][i] < r)
            for y in ind:
                p = np.random.rand()
                if p > 0.5:
                    Y[y,t+1] = 'I'

#%% SI Model with delay to infection

"""
This section of the code calculates the transition from S to I with a delay being considered in becoming infected since you can be infected and not infectious.

A matrix is calculated to determine what category each badger is in. This is done for each time point. Then, if the badger is within a certain distance (r)
for more than 20 times than the badger will make the transition from S to I.

S is susceptible and I is infectious.

"""

Y = np.zeros((n, nIter), dtype = str)

Y[:,0] = np.full((1, n), 'S', dtype=str)
Y[0,0] = 'I' #First badger is infected

C = np.zeros((n,1))

r = 0.25 #Distance to catch

for t in range (0,nIter-1):
    Y[:,t+1] = Y[:,t]
    for i in range (0,n):
        if Y[i,t] == 'I':
            ind = np.where((0 < D[t][0][i]) & (D[t][0][i] < r)) #Doesn't count itself
            for y in ind:
                C[y] +=1
        if C[i] > 20:
            Y[i,t+1] = 'I'
        

#%% Animation

"""
Then to view the transitions, we create an animation.
"""

fig = plt.figure(1)
fig.clf()
ims = []
for t in range(5000):
    c = ['#0000ff' if Y[i, t] == 'S' else '#ff0000' for i in range(n)]
    im = plt.scatter(X[::2,t], X[1::2,t], s=20, c=c)
    
    im.axes.set_xlim(-2.1, 2.1)
    im.axes.set_ylim(-2.1, 2.1)
    ims.append([im])
    
plt.axvline(x = 0, ymin = -3, ymax = 3, ls = '--', label = 'Border Line') #The boarder line (x=0)
# plt.axhline(y = 0, xmin = -4, xmax = 3, ls = '--') #The boarder line (y=0)

#Double Well
plt.plot(1,0,'go', label = 'P.O.A., i.e. Sett')
plt.plot(-1,0,'go')

#Quadruple Well
# plt.plot(1,-1,'go', label = 'P.O.A., i.e. Sett')
# plt.plot(-1,-1,'go')
# plt.plot(-1,1,'go')
# plt.plot(1,1,'go')

plt.legend(loc='lower left')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

anim = animation.ArtistAnimation(fig, ims, interval = 100, blit = True)
anim.save('badgers_infection.mp4')

#%% Animation to differentiate male and female

"""
If considering a more ecology approach within the code, the animation helps to differentiate between chosen females and males.
"""

fig = plt.figure(1)
fig.clf()
# fig.set_size_inches(5, 4, True)
ims = []
for t in range(1000):
    c1 = ['#0000ff' if Y[i, t] == 'S' else '#ff0000' for i in range(0,f)]
    c2 = ['#0000ff' if Y[i, t] == 'S' else '#ff0000' for i in range(f,f+m)]
    im1 = plt.scatter(X[0:2*f:2, t], X[1:2*f+1:2, t], s=30, c=c1 , marker = (5,2))
    im = plt.scatter(X[2*f:d:2, t], X[2*f+1:d:2, t], s=20 ,c=c2)
    im1.axes.set_xlim(-4, 3)
    im1.axes.set_ylim(-2.5, 2.5)
    ims.append([im1, im])

plt.axvline(x = 0, ymin = -3, ymax = 3, ls = '--', label = 'Border Line') #The boarder line (x=0)
# plt.axhline(y = 0, xmin = -4, xmax = 3, ls = '--') #The boarder line (y=0)

#Double Well
plt.plot(1,0,'go', label = 'P.O.A., i.e. Sett')
plt.plot(-1,0,'go')

#Quadruple Well
# plt.plot(1,-1,'go', label = 'P.O.A., i.e. Sett')
# plt.plot(-1,-1,'go')
# plt.plot(-1,1,'go')
# plt.plot(1,1,'go')

plt.legend(loc='lower left')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

anim = animation.ArtistAnimation(fig, ims, interval = 100, blit = True)
anim.save(f'badgersmf_{date}_infection.mp4')