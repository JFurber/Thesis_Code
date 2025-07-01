#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code investigates multiple badgers in a double-well energy potential with some ecology meanings.

Used in the chapter  'Modelling Strategies', under the section ' Ecology Model'

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Population(object):
    def __init__(self, S, epsilon_m, epsilon_f, alpha, beta):
        self.S = S # behaviour matrix
        self.n = S.shape[0]
        self.epsilon_m = epsilon_m # male noise value
        self.epsilon_f = epsilon_f # female noise value
        self.alpha = alpha # weight of attraction potential
        self.beta = beta #weight of repulsion potential

    def Attraction(self, x, i, j): 
        return -2*(x[2*i:2*i+2] - x[2*j:2*j+2])
    
    def Repulsion(self, x, i, j):
        return np.array([ 2*(x[2*i] - x[2*j])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2, 2*(x[2*i+1] - x[2*j+1])/((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2)**2])        

    def b(self, x):
        y = np.zeros_like(x) # Returns an array of given shape and type as the given array (i.e. x), with zeros
        for i in range(0, self.n):
            y[2*i:2*i+2] =  np.array([-4*x[2*i]**3 + 4*x[2*i] , -2*x[2*i+1]]) #Double well potential
            for j in range(0, self.n):
                if self.S[i][j] == 1:
                    y[2*i:2*i+2] += self.alpha * self.Attraction(x, i, j)
                elif self.S[i][j] == -1:
                    y[2*i:2*i+2] += self.beta * self.Repulsion(x, i, j)
        return y
        
    def sigma(self, x): #This is set for 10 badgers, with the first 4 being females and the next 6 being males.
        return np.diag([self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_f, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m, self.epsilon_m])

#%%
n = 10 # How many badgers are in the system

# %% Generation of S (the behaviour matrix)

"""
The general working of the ecology-focused behaviour matrix goes as follows:
The first f/2 females will be assigned to the first territory and the second half are part of the second territory. 
Then, with a similar approach to the males, the first m/2 are part of the first territory and the remainder are part of the second.
Then, it seems reasonable to assume that all females are attracted to all males, for mating purposes, regardless of social group.
In return, all males are attracted to all females. Then, if the badger shares the same territory, then
they will be attracted to them (since they are a social group) and all the others shall be repulsed.
"""

S1=np.ones((2,2))-np.eye(2)
S2=-1*np.ones((2,2))
S3=np.ones((3,3))-np.eye(3)
S4=-1*np.ones((3,3))
S5=np.ones((4,6))

S6=np.hstack((S1,S2))
S7=np.hstack((S2,S1))

S8=np.vstack((S6,S7,S5.T))

S9=np.hstack((S3,S4))
S10=np.hstack((S4,S3))

S11=np.vstack((S5,S9,S10))

S=np.hstack((S8,S11))

del S1, S2 , S3 , S4 , S5 , S6 , S7 , S8 , S9 , S10 , S11

# %%
P = Population(S, epsilon_m = 0.9, epsilon_f = 0.6 ,alpha = 0.01, beta = 0.01)
#%%
h = 1e-3 # step size
d = 2*n # State space
nIter = 5000 # number of iterations

X = np.zeros((d, nIter))
X[:,0] = 4*np.random.rand(d)-2 #the initial positions are randomised

#%% Euler-Maruyama and Plot
for t in range(1, nIter):
    X[:,t] = X[:,t-1] + P.b(X[:,t-1])*h + P.sigma(X[:,t-1])@np.random.randn(d)*np.sqrt(h)

plt.figure()
plt.clf()
for i in range(n):
    plt.plot(X[2*i,:], X[2*i+1,:])

plt.plot(X[0:8:2,0], X[1:9:2,0],'rx') # initial position of females 
plt.plot(X[8:20:2,0], X[9:20:2,0],'bx') # initial position of males
    
plt.axvline(x = 0, ymin = -3, ymax = 3, ls = '--', label = 'Border Line') #The boarder line

#Double Well
plt.plot(1,0,'go', label = 'P.O.A., i.e. Sett')
plt.plot(-1,0,'go')

plt.legend(loc='upper right')
plt.xlabel('Longitude')
plt.ylabel('Latitude')


#save the trajectory to plot in MATLAB (Check working directory before saving)
# np.savetxt('data.csv',X,delimiter=',')

#%% Animation

"""
The following code creates an animation for the trajectories. This can also be done on Matlab.
"""

ims = []
fig = plt.figure(1)
fig.clf()
for t in range(1000):
    im1, = plt.plot(X[0:8:2, t], X[1:9:2, t],'m.',animated=True)
    im, = plt.plot(X[8:20:2, t], X[9:20:2, t],'b.',animated=True)
    im.axes.set_xlim(-2.1, 2.1)
    im.axes.set_ylim(-2.1, 2.1)
    ims.append([im1,im])

anim = animation.ArtistAnimation(fig, ims, interval = 100, blit = True)
anim.save('badgers_ecology.mp4')

