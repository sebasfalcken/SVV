# -*- coding: utf-8 -*-

import numpy as np

C_a     = 0.484         #[m]
l_a     = 1.691         #[m]
x_1     = 0.149         #[m]
x_2     = 0.554         #[m]
x_3     = 1.541         #[m]
x_a     = 0.272         #[m]
h_a     = 0.173         #[m]
d_1     = 0.00681       #[m]
d_3     = 0.0203        #[m]
theta   = 26*np.pi/180  #[rad]
P       = 37900         #[N]
E       = 73100000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 
G       = 28000000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 

#Parameters that should be taken from others code:
z_tilde = -0.09056059254067199      #[m] 
I_zz    = 5.815938957599147e-06     #[m^4] 
I_yy    = 4.4539232840124055e-05    #[m^4]
J       = I_zz+I_yy                 #[m^4]
int4    = 2                         #[N*m^3]


#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5   and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12                   13    14        
#  0     1      2    3      4       5       6    7    8    9    10   11                   12    13

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

def dist(x):
    d   = np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])#14 elements
    tf  = d>=0
    return d*tf

def T(x):
    T_m = np.array([0, h_a/2*np.cos(theta)-(h_a/2-z_tilde)*np.sin(theta), 0, 0, z_tilde, z_tilde, z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)-(h_a/2-z_tilde)*np.sin(theta)), -1])
    
    
M_y     = np.array([-l_a+x_1, -np.cos(theta)*(l_a-x_2+x_a/2), -l_a+x_2, -l_a+x_3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
M_y_r   = -P*np.cos(theta)*(l_a-x_2-x_a/2)

M_z     = np.array([0, np.sin(theta)*(l_a-x_2+x_a/2), 0, 0, l_a-x_1, l_a-x_2, l_a-x_3, 0, 0, 0, 0, 0, 0, 0])
M_z_r   = P*np.sin(theta)*(l_a-x_2-x_a/2) 
#I still need to add the result of the integration

S_y     = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
S_y_r   = P*np.sin(theta)
#I still need to add the result of the integration

S_z     = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
S_z_r   = -P*np.cos(theta)

T       = np.array([0, h_a/2*np.cos(theta)-(h_a/2-z_tilde)*np.sin(theta), 0, 0, z_tilde, z_tilde, z_tilde, 0, 0, 0, 0, 0, 0, 0])

#So here I am puting in vector form the deflection and rotation equations of y, z and theta.
#There will be two additional values compared to the ones that I already put P and the integral. 
#These two values should be later separated and added to the right hand side of the equation

def v_def(x): 
    v       = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, 1/6, 1])
    #The last value should have Zeyad's function of integration in x
    v_2     = -1/E/I_zz*v*dist(x)**3 + np.array([0,0,0,0,0,0,0,x,1,0,0,0,0,0])
    v_2[12] = v_2[12] + v_2[13]
    return v_2[:-1] #Last value is the right hand side of the equation

def w_def(x):
    w       = np.array([-1/6, -np.cos(theta)/6, -1/6, -1/6, 0, 0, 0, 0, 0, 0, 0, 0, -1/6,0])
    return -1/E/I_yy*w*dist(x)**3 + np.array([0,0,0,0,0,0,0,0,0,x,1,0,0,0])

def th_rot(x):
    th      = np.array([0, h_a/2*np.cos(theta)-(h_a/2-z_tilde)*np.sin(theta), 0, 0, z_tilde, z_tilde, z_tilde, 0, 0, 0, 0, 0, -(-h_a/2*np.cos(theta)+(h_a/2-z_tilde)*np.sin(theta)), -1])
    th_2    = 1/G/J*th*dist(x) + np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0]) #vector of 14 elements
    th_2[12]= th_2[12]+th_2[13]
    return  th_2[:-1] #Last value is the right hand side of the equation

#Now, I will start with the boundary conditions

BC_vx1      = v_def(x_1)-z_tilde*th_rot(x_1)
BC_vx1_r    = BC_vx1[-1]
BC_vx1      = BC_vx1[:12]

BC_vx2      = v_def(x_2)-z_tilde*th_rot(x_2)
BC_vx2_r    = BC_vx2[-1]-d_1
BC_vx2      = BC_vx2[:12]

BC_vx3      = v_def(x_3)-z_tilde*th_rot(x_3)
BC_vx3_r    = BC_vx3[-1]-d_1+d_3
BC_vx3      = BC_vx3[:12]

BC_wx1      = w_def(x_1)
BC_wx1_r    = BC_wx1[-1]
BC_wx1      = BC_wx1[:12]

BC_wx2      = w_def(x_2)
BC_wx2_r    = BC_wx2[-1]
BC_wx2      = BC_wx2[:12]

BC_wx3      = w_def(x_3)
BC_wx3_r    = BC_wx3[-1]
BC_wx3      = BC_wx3[:12]

BC_wxac1    = w_def(x_2-x_a/2)
BC_wxac1_r  = BC_wxac1[-1]
BC_wxac1    = BC_wx3[:12]

#Now, the equations will be put in a matrix and will be solved :) :) !!!

A = np.array([M_y, M_z,S_y,S_z,T])

