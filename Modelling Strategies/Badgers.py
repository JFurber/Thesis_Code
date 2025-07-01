#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code investigates the multiple badgers in a double-well energy potential.

Used in the chapter  'Modelling Strategies', under the section 'Badger Model'

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Population(object):
    def __init__(self, S, sigma):
        self.S = S #the behaviour matrix
        self.n = S.shape[0]

    def Attraction(self, x, i, j): 
        return -2*(x[2*i:2*i+2] - x[2*j:2*j+2])
    
    def Repulsion(self, x, i, j):
        return np.array([ 2*(x[2*i] - x[2*j])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2, 2*(x[2*i+1] - x[2*j+1])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2])
        

    def b(self,x):
        y = np.zeros_like(x) # Returns an array of given shape and type as the given array (i.e. x), with zeros
        for i in range(0, self.n):
            y[2*i:2*i+2] =  np.array([-4*x[2*i]**3 + 4*x[2*i] , -2*x[2*i+1]]) #Double well potential
            for j in range(0, self.n):
                if self.S[i][j] == 1:
                    y[2*i:2*i+2] += 0.01 * self.Attraction(x, i, j)
                elif self.S[i][j] == -1:
                    y[2*i:2*i+2] += 0.01 * self.Repulsion(x, i, j)
        return y
        
    def sigma(self):
        return sigma

#%% Set up - Behaviour and Noise matrices
n = 5 # How many badgers are in the system

S=2*np.random.randint(0,2,size =(n,n))-1 #randomised attraction/repulsion matrix
np.fill_diagonal(S, 0)
print(S)

#chosen values for S
# S = np.array([[ 0,  1, -1, -1, -1], [ 1,  0,  1,  1, -1], [-1,  1,  0, -1,  1], [-1, -1,  1,  0, -1], [ 1, -1, -1, -1,  0]])


N = 0.9 #value of the noise

sigma = N * np.eye(2*n, 2*n) #noise matrix

P = Population(S,sigma)

h = 1e-3 # step size
d = 2*n # State space
nIter = 50000 # number of iterations
X = np.zeros((d, nIter))
X[:,0] = 4*np.random.rand(d)-2 #Initial values (random)


#%% Euler-Maruyama and Plot
for t in range(1, nIter):
    X[:,t] = X[:,t-1] + P.b(X[:,t-1])*h + P.sigma()@np.random.randn(d)*np.sqrt(h)

plt.figure(1)
plt.clf()
for i in range(n):
    plt.plot(X[2*i,:5000], X[2*i+1,:5000])
    plt.plot(X[2*i,0], X[2*i+1,0],'x')
    
#Save to plot in MATLAB (check working directory)
# np.savetxt('data.csv',X.T,delimiter=',')


#%% Animation
"""
The following code creates an animation for the trajectories. This can also be done on Matlab.
"""

ims = []
fig = plt.figure(1)
fig.clf()
for t in range(50):
    im, = plt.plot(X[::2, t], X[1::2, t],'b.',animated=True)
    im.axes.set_xlim(-2.1, 2.1)
    im.axes.set_ylim(-2.1, 2.1)
    ims.append([im])

anim = animation.ArtistAnimation(fig, ims, interval = 100, blit = True)
anim.save('badgers.mp4')
