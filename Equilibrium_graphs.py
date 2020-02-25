# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:05:28 2020

@author: Sebastian Falcken
"""

import numpy as np
import Equilibrium
import matplotlib.pyplot as plt

N=10000
x_i         = np.linspace(0,Equilibrium.l_a,N)

My          = []
for i in range(N):
    My      = np.append(My, Equilibrium.M_y(x_i[i]))

Mz          = []
for i in range(N):
    Mz      = np.append(Mz, Equilibrium.M_z(x_i[i]))

Tr          = []
for i in range(N):
    Tr      = np.append(Tr, Equilibrium.T_r(x_i[i]))

Sy          = []
for i in range(N):
    Sy      = np.append(Sy, Equilibrium.S_y(x_i[i]))

Sz          = []
for i in range(N):
    Sz      = np.append(Sz, Equilibrium.S_z(x_i[i]))

v           = []
for i in range(N):
    v       = np.append(v, Equilibrium.v_def(x_i[i]))

w           = []
for i in range(N):
    w       = np.append(w, Equilibrium.w_def(x_i[i]))

theta       = []
for i in range(N):
    theta   = np.append(theta, Equilibrium.th_rot(x_i[i]))


plt.figure()
plt.plot(x_i,My)
plt.show()

plt.figure()
plt.plot(x_i,Mz)
plt.show()

plt.figure()
plt.plot(x_i,Tr)
plt.show()

plt.figure()
plt.plot(x_i,Sy)
plt.show()

plt.figure()
plt.plot(x_i,Sz)
plt.show()

plt.figure()
plt.plot(x_i,v)
plt.show()

plt.figure()
plt.plot(x_i,w)
plt.show()

plt.figure()
plt.plot(x_i,theta)
plt.show()
