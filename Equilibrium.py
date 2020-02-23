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
E       = 73100000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 consulted
G       = 28000000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 consulted

#Parameters that should be taken from others code:
z_tilde = -0.09056059254067199      #[m] 
I_zz    = 5.815938957599147e-06     #[m^4] 
I_yy    = 4.4539232840124055e-05    #[m^4]
J       = I_zz+I_yy                 #[m^4]
int4    = 2                         #[N*m^3]


#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5  b_0,y b_0,2    and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12   13    14                      15  16

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

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

T       = np.array([0, -h_a/2*np.cos(theta)+(h_a/2-z_tilde)*np.sin(theta), 0, 0, z_tilde, z_tilde, z_tilde, 0, 0, 0, 0, 0, 0, 0])
#Still need to get results from Mitch

def v_def(x): #So here I am creating a vector which will have two additional values compared to the ones that I already put:
#P and the integral. These two values should be later separated and added to the right hand side of the equation
    dist    = np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])
    pos     = dist>=0
    v       = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, 0, 0, 1/6, 1])
    #The last value of v1 should have Zeyad's function of integration in x
    return -1/E/I_zz*pos*v*dist**3 + np.array([0,0,0,0,0,0,0,x,1,0,0,0,0,0,0,0])    

def w_def(x):
    dist    = np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, 0, 0, x-x_2-x_a/2])
    pos     = dist>=0
    w       = np.array([-1/6, -np.cos(theta)/6, -1/6, -1/6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/6])
    return -1/E/I_yy*pos*w*dist**3 + np.array([0,0,0,0,0,0,0,0,0,x,1,0,0,0,0])

def th_rot(x):
    dist    = np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])
    pos     = dist>=0
    th      = np.array([0, -h_a/2*np.cos(theta)+(h_a/2-z_tilde)*np.sin(theta), 0, 0, z_tilde, z_tilde, z_tilde, 0, 0, 0, 0, 0, 0, 0, -(h_a/2*np.cos(theta)-(h_a/2-z_tilde)*np.sin(theta)), -1])
    return 1/G/J*pos*th*dist + np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])

print(th_rot(x_3))
