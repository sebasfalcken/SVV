# -*- coding: utf-8 -*-
import numpy as np
import Properties_validation
import Aerodynamic_Load_Validation
#%% CASE 1: BENDING OF THE AILERON WITHOUT ANY LOADS APPLIED

C_a     = 0.605         #[m]
l_a     = 2.661         #[m]
x_1     = 0.172         #[m]
x_2     = 1.211         #[m]
x_3     = 2.591         #[m]
x_a     = 0.35          #[m]
h_a     = 0.205         #[m]
d_1     = 0.01154       #[m]
d_3     = 0.0184        #[m]
theta   = 28*np.pi/180  #[rad]
P       = 97400         #[N]
E       = 73100000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 
G       = 28000000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 

#Parameters that should be taken from others code:

I_zz, I_yy  = Properties_validation.test.MOI()
J           = Properties_validation.test.torsional_stiffness()  #[m^4]
z_tilde     = Properties_validation.test.Shear_center()[0]        #[m] 

#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5   and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12                   13    14        
#  0     1      2    3      4       5       6    7    8    9    10   11                   12    13

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

def dist(x): #Distance matrix for Macaulay functions
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0])

def dist1(x): #Distance matrix when you add the P value and the integrator
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])

M_y_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])*dist(l_a) 
M_y_t_r = -P*np.cos(theta)*(l_a-x_2-x_a/2)

M_z_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])*dist(l_a)
M_z_t_r = P*np.sin(theta)*(l_a-x_2-x_a/2)

S_z_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])
S_z_t_r = -P*np.cos(theta)

S_y_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])
S_y_t_r = P*np.sin(theta)

T_t     = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0])
T_t_r   = (h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P

def v_def_t(x):
    v_a     = -np.array([0, np.sin(theta)/6,0,0,1/6,1/6,1/6,0,0,0,0,0,P*np.cos(theta)/6, 0])/E/I_zz*(dist1(x)>=0)*dist1(x)**3
    v_a[7]  = x 
    v_a[8]  = 1
    return (v_a[:12], np.sum(v_a[-2:]))

def w_def_t(x):
    w_a     = -np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,-P*np.sin(theta)/6,0])/E/I_yy*(dist1(x)>=0)*dist1(x)**3
    w_a[9]  = x 
    w_a[10] = 1
    return (w_a[:12], np.sum(w_a[-2:]))

def th_rot_t(x):
    th_a    = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), 0])/G/J*(dist1(x)>=0)*dist1(x)
    th_a[11]= 1
    return (th_a[:12], np.sum(th_a[-2:]))

BC1     = v_def_t(x_1)[0]+th_rot_t(x_1)[0]*(z_tilde+h_a/2)
BC_1_r  = v_def_t(x_1)[1]+th_rot_t(x_1)[1]*(z_tilde+h_a/2) + d_1*np.cos(theta)

BC2     = v_def_t(x_2)[0]+th_rot_t(x_2)[0]*(z_tilde+h_a/2)
BC_2_r  = v_def_t(x_2)[1]+th_rot_t(x_2)[1]*(z_tilde+h_a/2)

BC3     = v_def_t(x_3)[0]+th_rot_t(x_3)[0]*(z_tilde+h_a/2)
BC_3_r  = v_def_t(x_3)[1]+th_rot_t(x_3)[1]*(z_tilde+h_a/2) +d_3*np.cos(theta)

BC4     = w_def_t(x_1)[0]
BC_4_r  = w_def_t(x_1)[1]-d_1*np.sin(theta)

BC5     = w_def_t(x_2)[0]
BC_5_r  = w_def_t(x_2)[1]

BC6     = w_def_t(x_3)[0]
BC_6_r  = w_def_t(x_3)[1]-d_3*np.sin(theta)

BC7     = w_def_t(x_2-x_a/2)[0]*np.cos(theta)+v_def_t(x_2-x_a/2)[0]*np.sin(theta)+th_rot_t(x_2-x_a/2)[0]*z_tilde*np.sin(theta)
BC_7_r  = w_def_t(x_2-x_a/2)[1]*np.cos(theta)+v_def_t(x_2-x_a/2)[1]*np.sin(theta)+th_rot_t(x_2-x_a/2)[1]*z_tilde*np.sin(theta)

A   = np.array([M_y_t, M_z_t, S_y_t, S_z_t, T_t, BC1, BC2, BC3, BC4, BC5, BC6, BC7])
b   = np.array([M_y_t_r, M_z_t_r, S_y_t_r, S_z_t_r, T_t_r, BC_1_r, BC_2_r, BC_3_r, BC_4_r, BC_5_r, BC_6_r, BC_7_r])
R   = np.linalg.solve(A,b)
R1  = np.append(R,[1,1])

def M_y(x):
    M1  = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0, P*np.cos(theta),0])
    return np.sum(M1*dist1(x)*(dist1(x)>=0)*R1)

def M_z(x):
    M2  = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0, -P*np.sin(theta),0])
    return np.sum(M2*dist1(x)*(dist1(x)>=0)*R1)

def T_r(x):
    T_f     = np.array([0, h_a/2*np.cos(theta)+z_tilde*np.sin(theta), 0, 0, h_a/2+z_tilde, h_a/2+z_tilde, h_a/2+z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, 0])
    return np.sum((dist1(x)>=0)*T_f*R1)

def S_y(x):
    S_ym    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -P*np.sin(theta), 0])
    return np.sum((dist1(x)>=0)*S_ym*R1)

def S_z(x):
    S_zm    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, np.cos(theta)*P, 0])
    return np.sum((dist1(x)>=0)*S_zm*R1)

def v_def(x):
    v   = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, -1/6*np.sin(theta)*P, 0])
    return np.sum(-1/E/I_zz*R1*v*(dist1(x)>=0)*dist1(x)**3)+R1[7]*x+R1[8]
    
def w_def(x):
    w = np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,P*np.sin(theta)/6,0])
    return np.sum(-1/E/I_yy*R1*w*(dist1(x)>=0)*dist1(x)**3)+R1[9]*x+R1[10]

def th_rot(x):
    th= np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, -P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), 0])
    return np.sum(th/G/J*(dist1(x)>=0)*dist1(x)*R1)+R1[11]
'''
#%% CASE 2: CLASSICAL CASE

C_a     = 0.605         #[m]
l_a     = 2.661         #[m]
x_1     = 0.172         #[m]
x_2     = 1.211         #[m]
x_3     = 2.591         #[m]
x_a     = 0.35          #[m]
h_a     = 0.205         #[m]
d_1     = 0.01154       #[m]
d_3     = 0.0184        #[m]
theta   = 28*np.pi/180  #[rad]
P       = 97400         #[N]
E       = 73100000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 
G       = 28000000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 

#Parameters that should be taken from others code:

I_zz, I_yy  = Properties_validation.test.MOI()
J           = Properties_validation.test.torsional_stiffness()  #[m^4]
z_tilde     = Properties_validation.test.Shear_center()[0]        #[m] 

#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5   and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12                   13    14        
#  0     1      2    3      4       5       6    7    8    9    10   11                   12    13

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

def dist(x): #Distance matrix for Macaulay functions
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0])

def dist1(x): #Distance matrix when you add the P value and the integrator
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])

M_y_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])*dist(l_a) 
M_y_t_r = -P*np.cos(theta)*(l_a-x_2-x_a/2)

M_z_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])*dist(l_a)
M_z_t_r = P*np.sin(theta)*(l_a-x_2-x_a/2)+Aerodynamic_Load_Verification.q2(l_a)*1000

S_z_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])
S_z_t_r = -P*np.cos(theta)

S_y_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])
S_y_t_r = P*np.sin(theta)+Aerodynamic_Load_Verification.q1(l_a)*1000

T_t     = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0])
T_t_r   = (h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P+Aerodynamic_Load_Verification.tau1(l_a)*1000

def v_def_t(x):
    v_a     = -np.array([0, np.sin(theta)/6,0,0,1/6,1/6,1/6,0,0,0,0,0,P*np.cos(theta)/6, Aerodynamic_Load_Verification.q4(x)*1000])/E/I_zz*(dist1(x)>=0)*dist1(x)**3
    v_a[7]  = x 
    v_a[8]  = 1
    return (v_a[:12], np.sum(v_a[-2:]))

def w_def_t(x):
    w_a     = -np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,-P*np.sin(theta)/6,0])/E/I_yy*(dist1(x)>=0)*dist1(x)**3
    w_a[9]  = x 
    w_a[10] = 1
    return (w_a[:12], np.sum(w_a[-2:]))

def th_rot_t(x):
    th_a    = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), Aerodynamic_Load_Verification.tau2(x)*1000])/G/J*(dist1(x)>=0)*dist1(x)
    th_a[11]= 1
    return (th_a[:12], np.sum(th_a[-2:]))

BC1     = v_def_t(x_1)[0]+th_rot_t(x_1)[0]*(z_tilde+h_a/2)
BC_1_r  = v_def_t(x_1)[1]+th_rot_t(x_1)[1]*(z_tilde+h_a/2) + d_1*np.cos(theta)

BC2     = v_def_t(x_2)[0]+th_rot_t(x_2)[0]*(z_tilde+h_a/2)
BC_2_r  = v_def_t(x_2)[1]+th_rot_t(x_2)[1]*(z_tilde+h_a/2)

BC3     = v_def_t(x_3)[0]+th_rot_t(x_3)[0]*(z_tilde+h_a/2)
BC_3_r  = v_def_t(x_3)[1]+th_rot_t(x_3)[1]*(z_tilde+h_a/2) +d_3*np.cos(theta)

BC4     = w_def_t(x_1)[0]
BC_4_r  = w_def_t(x_1)[1]-d_1*np.sin(theta)

BC5     = w_def_t(x_2)[0]
BC_5_r  = w_def_t(x_2)[1]

BC6     = w_def_t(x_3)[0]
BC_6_r  = w_def_t(x_3)[1]-d_3*np.sin(theta)

BC7     = w_def_t(x_2-x_a/2)[0]*np.cos(theta)+v_def_t(x_2-x_a/2)[0]*np.sin(theta)+th_rot_t(x_2-x_a/2)[0]*z_tilde*np.sin(theta)
BC_7_r  = w_def_t(x_2-x_a/2)[1]*np.cos(theta)+v_def_t(x_2-x_a/2)[1]*np.sin(theta)+th_rot_t(x_2-x_a/2)[1]*z_tilde*np.sin(theta)

A   = np.array([M_y_t, M_z_t, S_y_t, S_z_t, T_t, BC1, BC2, BC3, BC4, BC5, BC6, BC7])
b   = np.array([M_y_t_r, M_z_t_r, S_y_t_r, S_z_t_r, T_t_r, BC_1_r, BC_2_r, BC_3_r, BC_4_r, BC_5_r, BC_6_r, BC_7_r])
R   = np.linalg.solve(A,b)
R1  = np.append(R,[1,1])

def M_y(x):
    M1  = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0, P*np.cos(theta),0])
    return np.sum(M1*dist1(x)*(dist1(x)>=0)*R1)

def M_z(x):
    M2  = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0, -P*np.sin(theta),-Aerodynamic_Load_Verification.q2(x)*1000])
    return np.sum(M2*dist1(x)*(dist1(x)>=0)*R1)

def T_r(x):
    T_f     = np.array([0, h_a/2*np.cos(theta)+z_tilde*np.sin(theta), 0, 0, h_a/2+z_tilde, h_a/2+z_tilde, h_a/2+z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, -Aerodynamic_Load_Verification.tau1(x)*1000])
    return np.sum((dist1(x)>=0)*T_f*R1)

def S_y(x):
    S_ym    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -P*np.sin(theta), -Aerodynamic_Load_Verification.q1(x)*1000])
    return np.sum((dist1(x)>=0)*S_ym*R1)

def S_z(x):
    S_zm    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, np.cos(theta)*P, 0])
    return np.sum((dist1(x)>=0)*S_zm*R1)

def v_def(x):
    v   = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, -1/6*np.sin(theta)*P, -Aerodynamic_Load_Verification.q4(x)*1000])
    return np.sum(-1/E/I_zz*R1*v*(dist1(x)>=0)*dist1(x)**3)+R1[7]*x+R1[8]
    
def w_def(x):
    w = np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,P*np.sin(theta)/6,0])
    return np.sum(-1/E/I_yy*R1*w*(dist1(x)>=0)*dist1(x)**3)+R1[9]*x+R1[10]

def th_rot(x):
    th= np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, -P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), Aerodynamic_Load_Verification.tau2(x)*1000])
    return np.sum(th/G/J*(dist1(x)>=0)*dist1(x)*R1)+R1[11]

#%% CASE 3: UNBENT AILERON

C_a     = 0.605         #[m]
l_a     = 2.661         #[m]
x_1     = 0.172         #[m]
x_2     = 1.211         #[m]
x_3     = 2.591         #[m]
x_a     = 0.35          #[m]
h_a     = 0.205         #[m]
d_1     = 0.01154       #[m]
d_3     = 0.0184        #[m]
theta   = 28*np.pi/180  #[rad]
P       = 97400         #[N]
E       = 73100000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 
G       = 28000000000   #[Pa] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3 

#Parameters that should be taken from others code:

I_zz, I_yy  = Properties_validation.test.MOI()
J           = Properties_validation.test.torsional_stiffness()  #[m^4]
z_tilde     = Properties_validation.test.Shear_center()[0]        #[m] 

#So, rows are equations and columns are variables, just like in linear algebra.

#The variables will go in this order: 
#R_1,z  R_I  R_2,z  R_3,z  R_1,y  R_2,y  R_3,y  c_1  c_2  c_3  c_4  c_5   and ocasionally  P  int
#  1     2      3    4      5       6       7    8    9    10   11   12                   13    14        
#  0     1      2    3      4       5       6    7    8    9    10   11                   12    13

#Please reffer to the latex for more information.

#Additional assumptions:
#- The change in theta from actuator I to actuator II is negligible for the P components calculation

def dist(x): #Distance matrix for Macaulay functions
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0])

def dist1(x): #Distance matrix when you add the P value and the integrator
    return np.array([x-x_1, x-x_2+x_a/2, x-x_2, x-x_3, x-x_1, x-x_2, x-x_3, 0, 0, 0, 0, 0, x-x_2-x_a/2, 1])

M_y_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])*dist(l_a) 
M_y_t_r = -P*np.cos(theta)*(l_a-x_2-x_a/2)

M_z_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])*dist(l_a)
M_z_t_r = P*np.sin(theta)*(l_a-x_2-x_a/2)+Aerodynamic_Load_Verification.q2(l_a)*1000

S_z_t   = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0])
S_z_t_r = -P*np.cos(theta)

S_y_t   = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0])
S_y_t_r = P*np.sin(theta)+Aerodynamic_Load_Verification.q1(l_a)*1000

T_t     = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0])
T_t_r   = (h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P+Aerodynamic_Load_Verification.tau1(l_a)*1000

def v_def_t(x):
    v_a     = -np.array([0, np.sin(theta)/6,0,0,1/6,1/6,1/6,0,0,0,0,0,P*np.cos(theta)/6, Aerodynamic_Load_Verification.q4(x)*1000])/E/I_zz*(dist1(x)>=0)*dist1(x)**3
    v_a[7]  = x 
    v_a[8]  = 1
    return (v_a[:12], np.sum(v_a[-2:]))

def w_def_t(x):
    w_a     = -np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,-P*np.sin(theta)/6,0])/E/I_yy*(dist1(x)>=0)*dist1(x)**3
    w_a[9]  = x 
    w_a[10] = 1
    return (w_a[:12], np.sum(w_a[-2:]))

def th_rot_t(x):
    th_a    = np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), Aerodynamic_Load_Verification.tau2(x)*1000])/G/J*(dist1(x)>=0)*dist1(x)
    th_a[11]= 1
    return (th_a[:12], np.sum(th_a[-2:]))

BC1     = v_def_t(x_1)[0]+th_rot_t(x_1)[0]*(z_tilde+h_a/2)
BC_1_r  = v_def_t(x_1)[1]+th_rot_t(x_1)[1]*(z_tilde+h_a/2) + d_1*np.cos(theta)

BC2     = v_def_t(x_2)[0]+th_rot_t(x_2)[0]*(z_tilde+h_a/2)
BC_2_r  = v_def_t(x_2)[1]+th_rot_t(x_2)[1]*(z_tilde+h_a/2)

BC3     = v_def_t(x_3)[0]+th_rot_t(x_3)[0]*(z_tilde+h_a/2)
BC_3_r  = v_def_t(x_3)[1]+th_rot_t(x_3)[1]*(z_tilde+h_a/2) +d_3*np.cos(theta)

BC4     = w_def_t(x_1)[0]
BC_4_r  = w_def_t(x_1)[1]-d_1*np.sin(theta)

BC5     = w_def_t(x_2)[0]
BC_5_r  = w_def_t(x_2)[1]

BC6     = w_def_t(x_3)[0]
BC_6_r  = w_def_t(x_3)[1]-d_3*np.sin(theta)

BC7     = w_def_t(x_2-x_a/2)[0]*np.cos(theta)+v_def_t(x_2-x_a/2)[0]*np.sin(theta)+th_rot_t(x_2-x_a/2)[0]*z_tilde*np.sin(theta)
BC_7_r  = w_def_t(x_2-x_a/2)[1]*np.cos(theta)+v_def_t(x_2-x_a/2)[1]*np.sin(theta)+th_rot_t(x_2-x_a/2)[1]*z_tilde*np.sin(theta)

A   = np.array([M_y_t, M_z_t, S_y_t, S_z_t, T_t, BC1, BC2, BC3, BC4, BC5, BC6, BC7])
b   = np.array([M_y_t_r, M_z_t_r, S_y_t_r, S_z_t_r, T_t_r, BC_1_r, BC_2_r, BC_3_r, BC_4_r, BC_5_r, BC_6_r, BC_7_r])
R   = np.linalg.solve(A,b)
R1  = np.append(R,[1,1])

def M_y(x):
    M1  = np.array([-1, -np.cos(theta),-1,-1,0,0,0,0,0,0,0,0, P*np.cos(theta),0])
    return np.sum(M1*dist1(x)*(dist1(x)>=0)*R1)

def M_z(x):
    M2  = np.array([0, np.sin(theta),0,0,1,1,1,0,0,0,0,0, -P*np.sin(theta),-Aerodynamic_Load_Verification.q2(x)*1000])
    return np.sum(M2*dist1(x)*(dist1(x)>=0)*R1)

def T_r(x):
    T_f     = np.array([0, h_a/2*np.cos(theta)+z_tilde*np.sin(theta), 0, 0, h_a/2+z_tilde, h_a/2+z_tilde, h_a/2+z_tilde, 0, 0, 0, 0, 0,-(h_a/2*np.cos(theta)+z_tilde*np.sin(theta))*P, -Aerodynamic_Load_Verification.tau1(x)*1000])
    return np.sum((dist1(x)>=0)*T_f*R1)

def S_y(x):
    S_ym    = np.array([0, np.sin(theta), 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -P*np.sin(theta), -Aerodynamic_Load_Verification.q1(x)*1000])
    return np.sum((dist1(x)>=0)*S_ym*R1)

def S_z(x):
    S_zm    = np.array([-1, -np.cos(theta), -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, np.cos(theta)*P, 0])
    return np.sum((dist1(x)>=0)*S_zm*R1)

def v_def(x):
    v   = np.array([0, np.sin(theta)/6, 0, 0, 1/6, 1/6, 1/6, 0, 0, 0, 0, 0, -1/6*np.sin(theta)*P, -Aerodynamic_Load_Verification.q4(x)*1000])
    return np.sum(-1/E/I_zz*R1*v*(dist1(x)>=0)*dist1(x)**3)+R1[7]*x+R1[8]
    
def w_def(x):
    w = np.array([-1/6, -np.cos(theta)/6,-1/6,-1/6,0,0,0,0,0,0,0,0,P*np.sin(theta)/6,0])
    return np.sum(-1/E/I_yy*R1*w*(dist1(x)>=0)*dist1(x)**3)+R1[9]*x+R1[10]

def th_rot(x):
    th= np.array([0,h_a/2*np.cos(theta)+z_tilde*np.sin(theta),0,0,z_tilde+h_a/2,z_tilde+h_a/2,z_tilde+h_a/2,0,0,0,0,0, -P*(h_a/2*np.cos(theta)+z_tilde*np.sin(theta)), Aerodynamic_Load_Verification.tau2(x)*1000])
    return np.sum(th/G/J*(dist1(x)>=0)*dist1(x)*R1)+R1[11]
'''