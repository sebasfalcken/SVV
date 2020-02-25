# -*- coding: utf-8 -*-

import numpy as np
import Properties
import Aerodynamic_Load

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

I_zz, I_yy  = Properties.test.MOI()
J           = Properties.test.torsional_stiffness()  #[m^4]
z_tilde     = Properties.test.Shear_center()[0]        #[m] 

#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5   and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12                   13    14        
#  0     1      2    3      4       5       6    7    8    9    10   11                   12    13

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

def dist(x):
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])#14 elements
    #This will return the Macaulay term for each of the variables according to an x value. 

def T_t(x):
    T_m     = np.array([0, -h_a/2*np.cos(theta)-z_tilde*np.sin(theta), 0, 0, -h_a/2-z_tilde, -h_a/2-z_tilde, -h_a/2-z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, -Aerodynamic_Load.tau1(x)*1000])
    T_f     = (dist(x)>=0)*T_m
    T_f[12] = T_f[12]+T_f[13]
    return T_f[:-1]

def M_y_t(x):
    M_ym    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -np.cos(theta)*P, 0])
    M_yf    = dist(x)*(dist(x)>=0)*M_ym
    M_yf[12]= M_yf[12]+M_yf[13]
    return M_yf[:-1]

def M_z_t(x):
    M_zm    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, P*np.sin(theta), Aerodynamic_Load.q2(x)*1000])
    M_zf    = dist(x)*(dist(x)>=0)*M_zm
    M_zf[12]= M_zf[12]+M_zf[13]
    return M_zf[:-1]

def S_y_t(x):
    S_ym    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, np.sin(theta)*P, Aerodynamic_Load.q1(x)*1000])
    S_yf    = (dist(x)>=0)*S_ym
    S_yf[12]= S_yf[12]+S_yf[13]
    return S_yf[:-1]

def S_z_t(x):
    S_zm    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -P*np.cos(theta), 0])
    S_zf    = (dist(x)>=0)*S_zm
    S_zf[12]= S_zf[12]+S_zf[13]
    return S_zf[:-1]

#So here I am puting in vector form the deflection and rotation equations of y, z and theta.
#There will be two additional values compared to the ones that I already put P and the integral. 
#These two values should be later separated and added to the right hand side of the equation

def v_def_t(x): 
    v       = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, 1/6*np.sin(theta)*P, Aerodynamic_Load.q4(x)*1000])
    v_2     = -1/E/I_zz*v*(dist(x)*(dist(x)>=0))**3 + np.array([0,0,0,0,0,0,0,x,1,0,0,0,0,0])
    v_2[12] = v_2[12] + v_2[13]
    return v_2[:-1] #Last value is the right hand side of the equation

def w_def_t(x):
    w       = np.array([-1/6, -np.cos(theta)/6, -1/6, -1/6, 0, 0, 0, 0, 0, 0, 0, 0, -1/6*P*np.cos(theta),0])
    return -1/E/I_yy*w*(dist(x)*(dist(x)>=0))**3 + np.array([0,0,0,0,0,0,0,0,0,x,1,0,0,0])

def th_rot_t(x):
    th      = np.array([0, -h_a/2*np.cos(theta)-z_tilde*np.sin(theta), 0, 0, -h_a/2-z_tilde, -h_a/2-z_tilde, -h_a/2-z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, -Aerodynamic_Load.tau2(x)*1000])
    th_2    = 1/G/J*th*dist(x)*(dist(x)>=0) + np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0]) #vector of 14 elements
    th_2[12]= th_2[12]+th_2[13]
    return  th_2[:-1] #Last value is the right hand side of the equation

#Now, I will start with the boundary conditions

BC_vx1      = v_def_t(x_1)-(h_a/2+z_tilde)*th_rot_t(x_1)
BC_vx1_r    = BC_vx1[-1]
BC_vx1      = BC_vx1[:12]

BC_vx2      = v_def_t(x_2)-(h_a/2+z_tilde)*th_rot_t(x_2)
BC_vx2_r    = BC_vx2[-1]-d_1
BC_vx2      = BC_vx2[:12]

BC_vx3      = v_def_t(x_3)-(h_a/2+z_tilde)*th_rot_t(x_3)
BC_vx3_r    = BC_vx3[-1]-d_1+d_3
BC_vx3      = BC_vx3[:12]

BC_wx1      = w_def_t(x_1)
BC_wx1_r    = BC_wx1[-1]
BC_wx1      = BC_wx1[:12]

BC_wx2      = w_def_t(x_2)
BC_wx2_r    = BC_wx2[-1]
BC_wx2      = BC_wx2[:12]

BC_wx3      = w_def_t(x_3)
BC_wx3_r    = BC_wx3[-1]
BC_wx3      = BC_wx3[:12]

BC_wxac1    = w_def_t(x_2-x_a/2)
BC_wxac1_r  = BC_wxac1[-1]
BC_wxac1    = BC_wxac1[:12]

#Now, the equations will be put in a matrix and will be solved :) :) !!!

A = np.array([M_y_t(l_a)[:-1], M_z_t(l_a)[:-1], S_y_t(l_a)[:-1], S_z_t(l_a)[:-1], T_t(l_a)[:-1], BC_vx1, BC_vx2, BC_vx3, BC_wx1, BC_wx2, BC_wx3, BC_wxac1])
bf = np.array([M_y_t(l_a)[-1], M_z_t(l_a)[-1], S_y_t(l_a)[-1], S_z_t(l_a)[-1], T_t(l_a)[-1], BC_vx1_r, BC_vx2_r, BC_vx3_r, BC_wx1_r, BC_wx2_r, BC_wx3_r, BC_wxac1_r])

R_f     = (np.linalg.solve(A,bf)) #Reaction forces
#print(R_f)
R_f_1   = np.append(R_f, np.array([1, 1]))

def th_rot(x):
    th      = np.array([0, -h_a/2*np.cos(theta)-z_tilde*np.sin(theta), 0, 0, -h_a/2-z_tilde, -h_a/2-z_tilde, -h_a/2-z_tilde, 0, 0, 0, 0, 0,(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, -Aerodynamic_Load.tau2(x)*1000])
    th_f    = 1/G/J*th*dist(x)*(dist(x)>=0)*R_f_1 #The integral and the P value are still present after this calculation is done
    th_f[11]= R_f[11]
    return np.sum(th_f)

def M_y(x):
    M_ym    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, np.cos(theta)*P, 0])
    return np.sum(dist(x)*(dist(x)>=0)*M_ym*R_f_1)

def M_z(x):
    M_zm    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -P*np.sin(theta), -Aerodynamic_Load.q2(x)*1000])
    return np.sum(dist(x)*(dist(x)>=0)*M_zm*R_f_1)

def S_y(x):
    S_ym    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -np.sin(theta)*P, -Aerodynamic_Load.q1(x)*1000])
    return np.sum((dist(x)>=0)*S_ym*R_f_1)

def S_z(x):
    S_zm    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, P*np.cos(theta), 0])
    return np.sum((dist(x)>=0)*S_zm*R_f_1)

def v_def(x):
    v   = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, -1/6*P*np.sin(theta), -Aerodynamic_Load.q4(x)*1000])
    return np.sum(-1/E/I_zz*R_f_1*v*(dist(x)>=0)*dist(x)**3)+R_f_1[7]*x+R_f_1[8]
    
def w_def(x):
    w = np.array([-1/6, -np.cos(theta)/6, -1/6, -1/6, 0, 0, 0, 0, 0, 0, 0, 0, 1/6*P*np.cos(theta),0])
    return np.sum(-1/E/I_yy*R_f_1*w*(dist(x)>=0)*dist(x)**3)+R_f_1[9]*x+R_f_1[10]
'''
print(M_y(l_a))
print(M_z(l_a))
print(M_z_t(l_a))
print((dist(l_a)>=0))
'''