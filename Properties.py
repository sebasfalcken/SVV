import numpy as np
import math as m

class Properties:
    
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
    
    def changevar(self, Ca=0.484,la=1.691,x1=0.149,x2=0.554,x3=1.541,xa=0.27,ha=0.173,tsk=1.1/1000,tsp=2.5/1000,tst=1.2/1000,hst=14./1000,wst=18./1000,nst=13,d1=0.0681,d3=0.0203,theta=m.radians(26),P=37.9*1000):
        Vallist = [Ca,la,x1,x2,x3,xa,ha,tsk,tsp,tst,hst,wst,nst]
        if min(Vallist)<0:
            names = ["Ca","la","x1","x2","x3","xa","ha","tsk","tsp","tst","hst","wst","nst"]
            Vallist = [Ca,la,x1,x2,x3,xa,ha,tsk,tsp,tst,hst,wst,nst]
            index = Vallist.index(min(Vallist))
            raise Exception(f"{names[index]} cannot be smaller than zero")
        if nst%2 == 0:
            raise Exception(f"Amount of stringers cannot be smaller than even")
            
        self.Ca =  Ca  # m
        self.la =  la  # m
        self.x1 =  x1 # m
        self.x2 =  x2  # m
        self.x3 =  x3 # m
        self.xa =  xa # m
        self.ha =  ha # m
        self.tsk = tsk  # m
        self.tsp = tsp  # m
        self.tst = tst  # m
        self.hst = hst   # m
        self.wst = wst   # m
        self.nst = nst  # -
        self.d1 =  d1# m
        self.d3 =  d3 # m
        self.theta = theta  # rad
        self.P = P  # N
        Val_list = np.array([["Ca",self.Ca],["la",self.la],["x1",self.x1],["x2",self.x2],["x3",self.x3],["xa",self.xa],["ha",self.ha],
                             ["tsk",self.tsk],["tsp",self.tsp],["tst",self.tst],["hst",self.hst],["wst",self.wst],["nst",self.nst],
                             ["d1",self.d1],["d3",self.d3],["theta",self.theta],["P",self.P]])
        return Val_list
    
    def St_plcmnt(self):
        
        circumfer = m.pi*(self.ha/2.) + 2*np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.)
        amountstring = (self.nst - 1)
        distancestringer = np.linspace(0,circumfer,int(amountstring+2))
        dist_stringers = [[0,0]]

        for i in range(1,int(self.nst/2)+1):
            if distancestringer[i] < 0.5*m.pi*(self.ha/2.):
                theta = 2*((distancestringer[i]/self.ha))
                zst, yst = (self.ha/2)*(1-m.cos(theta)), (self.ha/2)*m.sin(theta)
                dist_stringers.append([-zst,yst])

            elif distancestringer[i] == 0.5*m.pi*(self.ha/2.):
                zst, yst = self.ha/2,self.ha/2
                dist_stringers.append([-zst,yst])

            else:
                dist_left = 0.5*m.pi*(self.ha/2.) - distancestringer[i]
                length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.)
                theta = m.atan((self.ha/2)/(self.Ca - self.ha/2.))
                zst, yst = self.Ca - (length_schuin-abs(dist_left))*m.cos(theta), (length_schuin-abs(dist_left))*m.sin(theta)
                dist_stringers.append([-zst,yst])

        stringers_reversed = dist_stringers[::-1]
        
        for j in range(len(stringers_reversed)-1):
            dist_stringers.append([stringers_reversed[j][0],-stringers_reversed[j][1]])

        return np.array(dist_stringers)
    
    def Centroid(self):
        dist_stringers = self.St_plcmnt()
        y_coord = 0
        z_loc = (self.ha/2 - 2*(self.ha/2)/(m.pi))*m.pi*(self.ha/2)*self.tsk + (self.ha/2)*self.ha*self.tsp + (((abs((self.ha/2)-self.Ca)/2))+self.ha/2)*2*np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.)*self.tsk
        area_st = (self.hst + self.wst)*self.tst
        for i in dist_stringers[:int(self.nst/2 + 1)]:
            z_loc += -2*i[0]*area_st

        z_area = m.pi*(self.ha/2)*self.tsk + self.ha*self.tsp + 2*np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.)*self.tsk + self.nst*area_st
        z_coord = -z_loc/z_area 
        
        return z_coord, y_coord
   
    def MOI(self):
        z_coord, y_coord = self.Centroid()
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        dist_stringers = self.St_plcmnt()
        
        theta = m.atan((self.ha/2)/(self.Ca - self.ha/2.))
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.)

        area_st = (self.hst + self.wst)*self.tst

        # IYY analysis

        Iyy_steiner_semi =  ((-1*(self.ha/2 - 2*(self.ha/2)/(m.pi)) - z_coord )**2) * m.pi*(self.ha/2)*self.tsk
        Iyy_spar = self.tsp*self.ha*((-self.ha/2 - z_coord)**2)
        Iyy_inclined = 2*((length_schuin**3)*self.tsk*(m.cos(theta)**2)/12  + length_schuin*self.tsk*(-1*((self.Ca - self.ha/2)/2 + self.ha/2)-z_coord)**2)
        I_semi = (m.pi*self.tsk*(self.ha/2)**3)/2

        #IZZ analysis

        Izz_spar = (self.tsp*(self.ha**3))/12
        Izz_inclined = 2*((length_schuin**3)*self.tsk*(m.sin(theta)**2)/12 + length_schuin*self.tsk*((self.ha/4)**2)) 

        # stringer analysis

        Izz_stringer = 0
        Iyy_stringer = 0

        for i in dist_stringers:
            Izz_stringer += area_st*(i[1])**2
            Iyy_stringer += area_st*(i[0]-z_coord)**2

        # total moment of inertia
        Izz_tot = Izz_spar + Izz_inclined + Izz_stringer + I_semi
        Iyy_tot = Iyy_spar + Iyy_inclined + Iyy_stringer + Iyy_steiner_semi + I_semi  
        return Izz_tot, Iyy_tot

    def total_area(self):
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        area1 = np.pi*(self.ha/2)*self.tsk + self.ha*self.tsp
        area2= 2*(length_schuin*self.tsk)
        stringers = self.nst*(self.hst*self.tst+self.wst*self.tst)
        return area1 + area2 + stringers
    
    def torsional_stiffness(self): 
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        area1 = 0.5*np.pi*(self.ha/2)**2
        area2 = (self.Ca-(self.ha/2))*(self.ha/2) 

        x1 = 1/(2*area1)*((1/self.tsk)*m.pi*(self.ha/2)+(1/self.tsp)*self.ha)
        x2 = -1/(2*area1)*(self.ha/self.tsp)
        x3 = -1/(2*area2)*(self.ha/self.tsp)
        x4 = 1/(2*area2)*((1/self.tsk)*2*length_schuin + self.ha/self.tsp)

        matrix = np.array([[2*area1, 2*area2, 0],[x1, x2, -1],[x3, x4,-1]])
        b = [1, 0, 0]

        x = np.linalg.solve(matrix, b)


        d_dz = (x[2])
        J = 1/(d_dz)
        return J
    
    def Shear_center(self):
        def summation(start,stop):
            if start < 0:
                return 0
            B_i = 0
            area_st = (self.hst + self.wst)*self.tst

            for i in dist_stringers[:int(self.nst/2+1)]:
                if abs(start) <= abs(i[0]) and abs(i[0]) <= abs(stop):
                    B_i += area_st*i[1]

            return B_i

        def integrateSin(N, lowerbound, upperbound):
            def f(x):
                return m.sin(x)
            number = 0
            number1 = 0

            for i in range(1, N + 1):
                number += f(lowerbound + (i-(1/2))*((upperbound-lowerbound)/N))
            number1 = ((upperbound-lowerbound)/N)*number

            return number1

        def integrate(N, lowerbound, upperbound):
            def f(x):
                return x
            number = 0
            number1 = 0

            for i in range(1, N + 1):
                number += f(lowerbound + (i-(1/2))*((upperbound-lowerbound)/N))
            number1 = ((upperbound-lowerbound)/N)*number

            return number1

        def shear(shearForce, Izz, thickness, y, integrate, lowerbound, upperbound, start, stop, qb_0):
            return (-shearForce/Izz)*(thickness*y*integrate(10000, lowerbound, upperbound) + summation(start,stop)) + qb_0
    
        Izz, Iyy = self.MOI()
        length_schuin = np.sqrt((self.ha/2.)**2 + (self.Ca - self.ha/2.)**2.) 
        dist_stringers = self.St_plcmnt()

        q1_shear = shear(1,Izz,self.tsk,(self.ha/2)**2,integrateSin,0,m.pi/2,0,self.ha/2,0)
        #q2_shear = shear(1,Izz,self.tsp,1,integrate,0,self.ha/2,-5,-5,0)
        q2_shear = (1/Izz)*self.tsk*(self.ha/2)
        q3_shear = shear(1,Izz,self.tsk,self.ha/2 - (self.ha/2)/length_schuin,integrate,0,length_schuin,self.ha/2,self.Ca,q1_shear+q2_shear)
        q4_shear = shear(1,Izz,self.tsk,(self.ha/2)/length_schuin,integrate,0,length_schuin,self.ha/2,self.Ca,q3_shear)
        q5_shear = (1/Izz)*self.tsk*(-self.ha/2)
        q6_shear = shear(1,Izz,self.tsk,(self.ha/2)**2,integrateSin,-m.pi/2,0,0,self.ha/2,q4_shear-q5_shear)

        x1 = (self.ha/2)*((m.pi/2)*2) + self.ha
        x2 = -1*(self.ha)
        x3 = -1*(self.ha)
        x4 = self.ha + 2*length_schuin

        b1 = (self.ha/2)*(q1_shear*(m.pi/2) + q6_shear*(m.pi/2))  + -1*q2_shear*(self.ha/2) + -1*q5_shear*self.ha/2
        b2 = q2_shear*self.ha/2 + q5_shear*self.ha/2 + q3_shear*length_schuin + q4_shear*length_schuin
        
        b = [-b1,-b2]
        
        matrix = np.array([[x1,x2],[x3,x4]])
        
        X = np.linalg.solve(matrix, b)

        q1_shear_tot = (q1_shear + X[0])*(m.pi*self.ha/2*0.5)*self.ha/2
        q3_shear_tot = (q3_shear + X[1] )*((self.Ca-self.ha/2)/2 + self.ha/2)*length_schuin
        q4_shear_tot = (q4_shear + X[1] )*((self.Ca-self.ha/2)/2 + self.ha/2)*length_schuin
        q6_shear_tot = (q6_shear + X[0])*(m.pi*self.ha/2*0.5)*self.ha/2

        Moment = q1_shear_tot + q3_shear_tot + q4_shear_tot + q6_shear_tot
        
        return (Moment-self.ha/2), 0


if __name__ == "__main__":
    test = Properties(1)
    z_coord, y_coord = test.Centroid()
    print(f"The z coord of the centroid is {z_coord}")
    print(f"The y coord of the centroid is {y_coord}")
    print()
    stringers = test.St_plcmnt()
    Izz,Iyy = test.MOI()
    print(f"The Izz is {Izz}")
    print(f"The Iyy is {Iyy}")
    print()
    total_area = test.total_area()
    print(f"Total area is {total_area}")
    print()
    J = test.torsional_stiffness()
    print(f"The Torsional stiffness is {J}")
    print()
    z_shear, y_shear = test.Shear_center()
    print(f"The z coord of the shear center is {z_shear}")
    print(f"The y coord of the shear center is {y_shear}")
    print()
    print(f"The array with sttingers are {stringers}")
