# -*- coding: utf-8 -*-

import numpy as np
import Equilibrium_Boeing
import matplotlib.pyplot as plt
import Validationfiles

#%%

My          = []
for i in range(len(Validationfiles.x_HL)):
    My      = np.append(My, Equilibrium_Boeing.M_y(Validationfiles.x_HL[i]))

Mz          = []
for i in range(len(Validationfiles.x_HL)):
    Mz      = np.append(Mz, Equilibrium_Boeing.M_z(Validationfiles.x_HL[i]))

Tr          = []
for i in range(len(Validationfiles.x_HL)):
    Tr      = np.append(Tr, Equilibrium_Boeing.T_r(Validationfiles.x_HL[i]))

Sy          = []
for i in range(len(Validationfiles.x_HL)):
    Sy      = np.append(Sy, Equilibrium_Boeing.S_y(Validationfiles.x_HL[i]))

Sz          = []
for i in range(len(Validationfiles.x_HL)):
    Sz      = np.append(Sz, Equilibrium_Boeing.S_z(Validationfiles.x_HL[i]))

v           = []
for i in range(len(Validationfiles.x_HL)):
    v       = np.append(v, Equilibrium_Boeing.v_def(Validationfiles.x_HL[i]))
v   = v*1000

w           = []
for i in range(len(Validationfiles.x_HL)):
    w       = np.append(w, Equilibrium_Boeing.w_def(Validationfiles.x_HL[i]))
w   = w*1000

theta1       = []
for i in range(len(Validationfiles.x_HL)):
    theta1   = np.append(theta1, Equilibrium_Boeing.th_rot(Validationfiles.x_HL[i]))

#%%

plt.figure()
plt.scatter(Validationfiles.x_HL,w, label="Numerical values")
plt.scatter(Validationfiles.x_HL,np.array(Validationfiles.defl_z_HL_bending),label="Validation values",color='red',s=10)
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$w(x) [mm]$')
plt.title('Deflection in z - Bent aileron without aerodynamic load applied')
plt.legend()
plt.show

plt.figure()
plt.scatter(Validationfiles.x_HL,v, label="Numerical values")
plt.scatter(Validationfiles.x_HL,np.array(Validationfiles.defl_y_HL_bending),label="Validation values",color='red',s=10)
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$v(x) [mm]$')
plt.title('Deflection in y - Bent aileron without aerodynamic load applied')
plt.legend()
plt.show
'''
plt.figure()
plt.plot(x_i,w, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_z_HL_jam_bent),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$w(x) [mm]$')
plt.title('Deflection in z - Bent aileron with jammed actuator and aerodynamic load')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,v, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_y_HL_jam_bent),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$v(x) [mm]$')
plt.title('Deflection in y - Bent aileron with jammed actuator and aerodynamic load')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,w, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_z_HL_jam_straight),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$w(x)[mm]$')
plt.title('Deflection in z - Unbent aileron with jammed actuator and aerodynamic load')
plt.legend()
plt.show

plt.figure()
plt.plot(x_i,v, label="Numerical values")
plt.plot(x_i,np.array(Validationfiles.defl_y_HL_jam_straight),label="Validation values")
plt.xlabel("Aileron span [m]")
plt.ylabel(r'$v(x) [mm]$')
plt.title('Deflection in y - Unbent aileron with jammed actuator and aerodynamic load')
plt.legend()
plt.show
'''