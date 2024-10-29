# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 12:02:38 2023

@author: Michael
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

n = np.linspace(0,100,101)
r_tube = (10**2)/(2*np.pi)  #Each square of the lattice is 10nm, so to get the circumference of 
                            #tube, you need 10*10, if it's 10 squares.
F = 2*np.pi*np.sqrt(n*(10**2)/np.pi)
F_tube = [2*2*np.pi*r_tube]*len(n[33:])   #This is only assuming if the cluster were to wrap
#The lines are distinct from the wrapping transition point.
#IT is only showing that once it hits a cluster size, compared to the round cluster
#the wrapped cluster would by definition, have a shorter interface,
#But these values are probably below the wrapping transition.
print(4*np.pi*(r_tube**2)/100)

r_tube2 = (10*13)/(2*np.pi) #Using 13 square
from_here2 = 54
F_tube2 = [2*2*np.pi*r_tube2]*len(n[from_here2:]) 
print(4*np.pi*(r_tube2**2)/100)
#Oh, the values are pretty close. 

plt.plot(n[from_here2:],F_tube2, linewidth = 4, c = 'black', linestyle = 'dashed', label = r"Wrapped ($d_{tube} = 4.1 _ d_{IRE1}$)")
plt.plot(n[33:],F_tube, linewidth = 4, c = 'r',  linestyle = 'dotted', label = r"Wrapped ($d_{tube} = 3.2 _ d_{IRE1}$)")
plt.plot(n,F,linewidth = 4, c= 'b', label = "Round")

plt.yticks([0,100,200,300], [0,10,20,30])
plt.ylabel('Interface length ($d_{IRE1}$)', fontsize = 18)
plt.xlabel('Cluster size (proteins)', fontsize = 18)
#plt.title('Length of Cluster Interface')
plt.legend(prop={'size': 12}, frameon = False, loc = 'upper left')

plt.savefig('ClusterInterface2_Scaled.svg')