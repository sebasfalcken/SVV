# -*- coding: utf-8 -*-

import numpy as np
import Equilibrium
import matplotlib.pyplot as plt
import main

#%%

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

theta1       = []
for i in range(N):
    theta1   = np.append(theta1, Equilibrium.th_rot(x_i[i]))

#%%

plt.figure()
plt.plot(x_i,My, label="Numerical values")
plt.plot(x_i,main.M_y,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$M_{y}(x) \ [N \cdot m]$')
plt.title('Bending moment about y')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,Sz, label="Numerical values")
plt.plot(x_i,main.S_z,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$S_{z}(x) \ [N]$')
plt.title('Shear force in z')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,w, label="Numerical values")
plt.plot(x_i,main.w,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$w(x)$')
plt.title('Deflection in z')
plt.legend()
plt.show


plt.figure()
plt.plot(x_i,Mz, label="Numerical values")
plt.plot(x_i,main.M_z,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$M_{z}(x) \ [N \cdot m]$')
plt.title('Bending moment about z')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,Sy, label="Numerical values")
plt.plot(x_i,main.S_y,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$S_{y}(x) \ [N]$')
plt.title('Shear force in y')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,v, label="Numerical values")
plt.plot(x_i,main.v,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$v(x)$')
plt.title('Deflection in y')
plt.legend()
plt.show


plt.figure()
plt.plot(x_i,Tr, label="Numerical values")
plt.plot(x_i,main.T_1,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$T(x) \ [N \cdot m]$')
plt.title('Torque')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,theta1, label="Numerical values")
plt.plot(x_i,main.th,label="Verification values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$\theta (x) \ [rad]$')
plt.title('Twist in x')
plt.legend()
plt.show