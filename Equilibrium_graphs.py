# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:05:28 2020

@author: Sebastian Falcken
"""

import numpy as np
import Equilibrium
import matplotlib.pyplot as plt

N=10000
x_i     = np.linspace(0,Equilibrium.l_a,N)

My      = []
for i in range(N):
    My  = np.append(My, Equilibrium.M_y(x_i[i]))

Mz      = []
for i in range(N):
    Mz  = np.append(Mz, Equilibrium.M_z(x_i[i]))

plt.figure()
plt.plot(x_i,My)
plt.show()

plt.figure()
plt.plot(x_i,Mz)
plt.show()