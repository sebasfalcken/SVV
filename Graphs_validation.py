# -*- coding: utf-8 -*-

import numpy as np
import Equilibrium_Boeing
import matplotlib.pyplot as plt
import Validationfiles

#%%

N=10000
x_i         = np.linspace(0,Equilibrium_Boeing.l_a,N)

My          = []
for i in range(N):
    My      = np.append(My, Equilibrium_Boeing.M_y(x_i[i]))

Mz          = []
for i in range(N):
    Mz      = np.append(Mz, Equilibrium_Boeing.M_z(x_i[i]))

Tr          = []
for i in range(N):
    Tr      = np.append(Tr, Equilibrium_Boeing.T_r(x_i[i]))

Sy          = []
for i in range(N):
    Sy      = np.append(Sy, Equilibrium_Boeing.S_y(x_i[i]))

Sz          = []
for i in range(N):
    Sz      = np.append(Sz, Equilibrium_Boeing.S_z(x_i[i]))

v           = []
for i in range(N):
    v       = np.append(v, Equilibrium_Boeing.v_def(x_i[i]))

w           = []
for i in range(N):
    w       = np.append(w, Equilibrium_Boeing.w_def(x_i[i]))

theta1       = []
for i in range(N):
    theta1   = np.append(theta1, Equilibrium_Boeing.th_rot(x_i[i]))

#%%

plt.figure()
plt.plot(x_i,w, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_z_HL_bending),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$w(x)$')
plt.title('Deflection in z')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,v, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_y_HL_bending),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$v(x) [m]$')
plt.title('Deflection in y - Aileron bending without aerodynamic load applied')
plt.legend()
plt.show
