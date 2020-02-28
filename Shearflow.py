import numpy as np
import Properties
import math as m
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


class ShearFlow:
    def __init__(self, plane):
        self.plane = plane
        self.Ca = 0.484  # m
        self.la = 1.691  # m
        self.x1 = 0.149  # m
        self.x2 = 0.554  # m
        self.x3 = 1.541  # m
        self.xa = 0.27  # m
        self.ha = 0.173  # m
        self.tsk = 1.1/1000  # m
        self.tsp = 2.5/1000  # m
        self.tst = 1.2/1000  # m
        self.hst = 14./1000   # m
        self.wst = 18./1000   # m
        self.nst = 13  # -
        self.d1 = 0.0681  # m
        self.d3 = 0.0203  # m
        self.theta = m.radians(26)  # rad
        self.P = 37.9*1000  # N
        
    def summation(self,start,stop):
        dist_stringers = Properties.Properties(1).St_plcmnt()
        if start < 0:
            return 0
        B_i = 0
        area_st = (self.hst + self.wst)*self.tst
        for i in dist_stringers[:7]:
            if abs(start) <= abs(i[0]) and abs(i[0]) <= abs(stop):
                B_i += area_st*i[1]
        return B_i

    def summationyy(self,start,stop):
        dist_stringers = Properties.Properties(1).St_plcmnt()
        z_coord, y_coord = Properties.Properties(1).Centroid()
        if start == stop:
            return 0
        B_i = 0
        area_st = (self.hst + self.wst)*self.tst
        for i in dist_stringers[:7]:
            if abs(start) <= abs(i[0]) and abs(i[0]) <= abs(stop):
                B_i += area_st*(i[0]-z_coord)
        return B_i

    def integrateSin(self, N, lowerbound, upperbound):
        def f(x):
            return m.sin(x)
        number = 0
        number1 = 0
    
        for i in range(1, N + 1):
            number += f(lowerbound + (i-(1/2))*((upperbound-lowerbound)/N))
        number1 = ((upperbound-lowerbound)/N)*number
    
        return number1

    def integrate(self, N, lowerbound, upperbound):
        def f(x):
            return x
        number = 0
        number1 = 0
    
        for i in range(1, N + 1):
            number += f(lowerbound + (i-(1/2))*((upperbound-lowerbound)/N))
        number1 = ((upperbound-lowerbound)/N)*number
    
        return number1

    def integrateCos(self, N, lowerbound, upperbound):
        def f(x):
            return m.cos(x)
        number = 0
        number1 = 0
    
        for i in range(1, N + 1):
            number += f(lowerbound + (i-(1/2))*((upperbound-lowerbound)/N))
        number1 = ((upperbound-lowerbound)/N)*number
    
        return number1

    def shear(self, shearForce, Izz, thickness, y, integrate, lowerbound, upperbound, start, stop, qb_0):
        return (-shearForce/Izz)*(thickness*y*integrate(10000, lowerbound, upperbound) + self.summation(start,stop)) + qb_0

    def shearyy(self,shearForce, Izz, thickness, y, integrate, lowerbound, upperbound, start, stop, qb_0):
        return (-shearForce/Izz)*(thickness*y*integrate(10000, lowerbound, upperbound) + self.summationyy(start,stop)) + qb_0
    
    def Print_shear(self, Sy, Sx, T):
        Izz, Iyy = Properties.Properties(1).MOI()
        z_coord, y_coord = Properties.Properties(1).Centroid()
        dist_stringers = Properties.Properties(1).St_plcmnt()
        area_st = (self.hst + self.wst)*self.tst

        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        area1 = 0.5*np.pi*(self.ha/2)**2
        area2 = (self.Ca-(self.ha/2))*(self.ha/2) 

        q1_shear = self.shear(Sy,Izz,self.tsk,(self.ha/2)**2,self.integrateSin,0,m.pi/2,0,self.ha/2,0)
        q2_shear = self.shear(Sy,Izz,self.tsp,1,self.integrate,0,self.ha/2,-5,-5,0)
        q3_shear = self.shear(Sy,Izz,self.tsk,self.ha/2 - (self.ha/2)/length_schuin,self.integrate,0,length_schuin,self.ha/2,self.Ca,q1_shear+q2_shear)
        q4_shear = self.shear(Sy,Izz,self.tsk,(self.ha/2)/length_schuin,self.integrate,0,length_schuin,self.ha/2,self.Ca,q3_shear)
        q5_shear = self.shear(Sy,Izz,self.tsp,1,self.integrate,0,-self.ha/2,-5,-5,0)
        q6_shear = self.shear(Sy,Izz,self.tsk,(self.ha/2)**2,self.integrateSin,-m.pi/2,0,0,self.ha/2,q4_shear-q5_shear)
        
        x1 = (self.ha)*((m.pi/2)*2)/self.tsk + self.ha/self.tsp
        x2 = -1*(self.ha)/self.tsp
        x3 = -1*(self.ha)/self.tsp
        x4 = self.ha/self.tsp + 2*length_schuin/self.tsk
        
        b1 = (self.ha/2)*(q1_shear*(m.pi/2) + q6_shear*(m.pi/2))/self.tsk  + q2_shear*(self.ha/2)/self.tsp + q5_shear*(self.ha/2)/self.tsp
        b2 = q2_shear*(self.ha/2)/self.tsp + q5_shear*(self.ha/2)/self.tsp + q3_shear*(length_schuin)/self.tsk + q4_shear*length_schuin/self.tsk
                
        b = [-b1,-b2]
        matrix = np.array([[x1,x2],[x3,x4]])
        
        X2 = np.linalg.solve(matrix, b)
        
        x1 = 1/(2*area1)*((1/self.tsk)*m.pi*(self.ha/2)+(1/self.tsp)*self.ha)
        x2 = -1/(2*area1)*(self.ha/self.tsp)
        x3 = -1/(2*area2)*(self.ha/self.tsp)
        x4 = 1/(2*area2)*((1/self.tsk)*2*length_schuin + self.ha/self.tsp)
        
        matrix = np.array([[2*area1, 2*area2, 0],[x1, x2, -1],[x3, x4,-1]])
        b = [-T, 0, 0]
        
        X = np.linalg.solve(matrix, b)
        
        #region one
        theta = np.linspace(0,m.pi/2,num = 100)
        q1 =[]
        
        for i in range(len(theta)):
            a = (self.ha/2)-(self.ha/2)*m.cos(theta[i])
            shear1 = 0
            if a >= abs(dist_stringers[1][0]):
                shear1 = (-Sy/Izz)*area_st*dist_stringers[1][1] + (-Sx/Iyy)*area_st*(dist_stringers[1][0]-z_coord)
                
            shear1_zz = self.shear(Sy,Izz,self.tsk,(self.ha/2)**2, self.integrateSin, 0, theta[i], -5, -5,0)
            shear1_yy = (-Sx/Iyy)*self.tsk*-1*(self.ha/2 + z_coord)*(self.ha/2)*theta[i]
            shear1_yy += self.shearyy(Sx,Iyy,self.tsk,(self.ha/2)**2,self.integrateCos,0,theta[i],-5,-5,0)
            q1.append((shear1_zz + shear1_yy + shear1 + X[0] + X2[0]))
        
        #region two
        y = np.linspace(0,self.ha/2.,num = 100)
        q2 =[]
        for i in range(len(y)):
            q2_shear  = self.shear(Sy,Izz,self.tsp,1,self.integrate,0,y[i],-5,-5,0)
            q2_shear_yy = (-Sx/Iyy)*self.tsp*(-self.ha/2-z_coord)*y[i]
            q2.append(q2_shear + q2_shear_yy + X[0] - X[1] + X2[1]-X2[0])
            
        #region three
        s = np.linspace(0,m.sqrt((self.Ca-self.ha/2.)**2+(self.ha/2.)**2),num = 100)
        q3 =[]
        for i in range(len(s)):
            q3_shear = self.shear(Sy,Izz,self.tsk,-1*(self.ha/2)/length_schuin,self.integrate,0,s[i],self.ha/2,self.ha/2+s[i],0)
            q3_shear += (-Sy/Izz)*self.tsk*(self.ha/2)*s[i]
            q3_shear_yy = (-Sx/Iyy)*(self.tsk*(-self.ha/2-z_coord))*s[i]
            q3_shear_yy += self.shearyy(Sx,Iyy,self.tsk,-(self.Ca - self.ha/2)/length_schuin,self.integrate,0,s[i],self.ha/2,self.ha/2+s[i],0)
            q3.append(q3_shear + q3_shear_yy + q1[-1] + q2[-1] + X[1])
        
        #region four
        s = np.linspace(m.sqrt((self.Ca-self.ha/2.)**2+(self.ha/2.)**2),0,num = 100)
        q4 = []
        theta = m.acos((self.Ca-self.ha/2)/length_schuin)
        counter = 7
        for j in range(len(s)):
            q4_shear = self.shear(Sy,Izz,self.tsk,(-self.ha/2)/length_schuin,self.integrate,0,s[j],-5,-5,0)
            q4_shear_yy = (-Sx/Iyy)*(self.tsk*(-self.Ca-z_coord))*s[j]
            q4_shear_yy += self.shearyy(Sx,Iyy,self.tsk,(self.Ca - self.ha/2)/length_schuin,self.integrate,0,s[j],-5,-5,0)
            for i in dist_stringers[6:counter]:
                if abs(i[0]) >= self.ha/2+s[j]*m.cos(theta):
                    q4_shear += (-Sy/Izz)*area_st*i[1] + (-Sx/Iyy)*area_st*(i[0]-z_coord)
                    counter += 1
               
            q4.append(q4_shear + q4_shear_yy  + q4_shear_yy + X[1])
        
        #region five
        y = np.linspace(0,-self.ha/2.,num = 100)
        q5 =[]
        for i in range(len(y)):
            q5_shear  = self.shear(Sy,Izz,self.tsp,1,self.integrate,0,y[i],-5,-5,0)
            q5_shear_yy = (-Sx/Iyy)*self.tsp*(-self.ha/2-z_coord)*y[i]
            q5.append(q5_shear + q5_shear_yy + X[1] - X[0] + X2[1] - X2[0])  
         
        
        #region six
        theta = np.linspace(-m.pi/2,0,num = 100)
        q6 = []
        shear6 = 0
        for i in range(len(theta)):
            a = (self.ha/2)-(self.ha/2)*m.cos(abs(theta[i]))
            if a <= abs(dist_stringers[12][0]):
                shear6 = (-Sy/Izz)*area_st*dist_stringers[12][1] + (-Sx/Iyy)*area_st*(dist_stringers[12][0] - z_coord)
            
            shear6_zz = self.shear(Sy,Izz,self.tsk,(self.ha/2)**2, self.integrateSin,  0, theta[i],-5, -5,0)
            shear6_yy = (-Sx/Iyy)*self.tsk*(-1*(self.ha/2+z_coord)*self.ha/2)*theta[i]
            shear6_yy += self.shearyy(Sx,Iyy,self.tsk,(self.ha/2)**2,self.integrateCos, 0,theta[i],-5,-5,0)
            q6.append((shear6_zz + shear6_yy + shear6  + q4[-1] - q5[-1] + X2[0]))
            
        theta = np.linspace(0,m.pi/2,num = 100)
        x = []
        y = []
        for i in theta:
            zst, yst = (self.ha/2)*(1-m.cos(i)), (self.ha/2)*m.sin(i)
            x.append(-zst)
            y.append(yst)
        x = np.array(x)
        y = np.array(y)
        
        q2_y = np.linspace(0,self.ha/2.,num = 100)
        q2_z = np.full((100,1),-self.ha/2)
        
        x = np.append(x,q2_z)
        y = np.append(y,q2_y)
        
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        s = np.linspace(0,m.sqrt((self.Ca-self.ha/2.)**2+(self.ha/2.)**2),num = 100)
        alpha = m.atan((self.ha/2)/(self.Ca - self.ha/2.))
        
        for i in s:
            x = np.append(x,-1*(self.Ca-(length_schuin - i)*m.cos(alpha)))
            y = np.append(y,((length_schuin - i)*m.sin(alpha)))
        for i in s[::-1]:
            x = np.append(x,-1*(self.Ca-(length_schuin - i)*m.cos(alpha)))
            y = np.append(y,(-1*(length_schuin - i)*m.sin(alpha)))
        
        q2_y = np.linspace(0,-self.ha/2.,num = 100)
        q2_z = np.full((100,1),-self.ha/2)
        
        x = np.append(x,q2_z)
        y = np.append(y,q2_y)
        
        for i in theta[::-1]:
            zst, yst = (self.ha/2)*(1-m.cos(i)), (self.ha/2)*m.sin(i)
            x = np.append(x,-zst)
            y = np.append(y,-yst)
        
               
        x = np.array(x)
        y = np.array(y)
        data = np.array(q1) # first derivative
        data = np.append(data,q2)
        data = np.append(data,q3)
        data = np.append(data,q4[::-1])
        data = np.append(data,q5)
        data = np.append(data,q6)
        return data, x, y
    
    def Plot(self, data, x, y):
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        fig, axis = plt.subplots(1, 1)
        
        norm = plt.Normalize(data.min(), data.max())
        line_col = LineCollection(segments, cmap='rainbow', norm=norm)
        # Set the values used for colormapping
        line_col.set_array(data)
        line_col.set_linewidth(2)
        line = axis.add_collection(line_col)
        fig.colorbar(line, ax=axis, label=r'q = $N/m $')
        
        # Use a boundary norm instead
        
        axis.set_xlim(0.05, -0.5)
        axis.set_ylim(-0.15, 0.15)
        axis.set_xlabel('z [m]')
        axis.set_title('Shear flow distribution')
        axis.set_ylabel('y [m]')
        plt.show()


if __name__ == "__main__":
    #Sy = -4985.0292997
    #Sz = -73359.55813566
    #T = -1743.7425962
    
    Sy = 56576.054162008    
    Sz = -39021.131846804834
    T = 1291.9302730065012
    
    #Sy = -71552.22303513
    #Sx = 202054.57756203
    #T = -1184.42632407
    
    a = ShearFlow(1)
    data, x, y = a.Print_shear(Sy, Sz, T)
    a.Plot(data, x, y)
