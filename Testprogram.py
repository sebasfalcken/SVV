import numpy as np
import Stiffness
import math as m
import unittest
import Properties 

aircraft = "CRJ700"
Ca = 0.484  # m
la = 1.691  # m
x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 0.272  # m
ha = 0.173  # m
tsk = 1.1/1000  # m
tsp = 2.5/1000  # m
tst = 1.2/1000  # m
hst = 14./1000   # m
wst = 18./1000   # m
nst = 13  # -
d1 = 0.00681  # m
d3 = 0.02030  # m
theta = m.radians(26)  # rad
P = 37.9*1000  # N

def SSD(lhs, rhs):
    return np.sum((lhs - rhs)**2)

class VerifyNum(unittest.TestCase):
    
    def __init__(self):
        self.num = Properties.Properties(1)
        self.num.changevar(Ca,la,x1,x2,x3,xa,ha,tsk,tsp,tst,hst,wst,nst,d1,d3,theta,P)
        self.nst = nst
        self.crosssection = Stiffness.Crosssection(nst,Ca,ha,tsk,tsp,tst,hst,wst)
        self.crosssection.compute_bending_properties()

        self.stringers = self.crosssection.stcoord            # array containing stringer coordinates
        self.tot_area = self.crosssection.totarea            # total cross-section area
        self.y = self.crosssection.yc                 # y-coordinate of the centroid
        self.z = self.crosssection.zc                 # z-coordinate of the centroid
        self.Iyy_corr = self.crosssection.Iyy                # moment of inertia about y-axis
        self.Izz_corr = self.crosssection.Izz 
        self.crosssection.compute_shearcenter()   # Run the calculations
        self.crosssection.compute_torsionalstiffness()
        
        self.ysc = self.crosssection.ysc                 # y-coordinate of the centroid
        self.zsc = self.crosssection.zsc                 # z-coordinate of the centroid
        self.J = self.crosssection.J 
       
    def test_stringer(self):
        
        for i in range(self.nst):
            assert(SSD(self.stringers[i][0], self.num.St_plcmnt()[i][0])< 0.000001)
            assert(SSD(self.stringers[i][1], self.num.St_plcmnt()[i][1])< 0.000001)
    
    def test_crossextionalArea(self):
        area = self.num.total_area()
        assert(SSD(self.tot_area,area)<0.000001)
        
    def test_centroid(self):
        z_coord,y_coord = self.num.Centroid()
        assert(SSD(self.y, y_coord)< 0.000001)
        assert(SSD(self.z, z_coord)< 0.000001)
        
    def test_momentOfinertia(self):
        Izz, Iyy = self.num.MOI()
        assert(SSD(self.Izz_corr, Izz)< 0.000001)
        assert(SSD(self.Iyy_corr, Iyy)< 0.000001)
    
    def test_torsional_stiffenss(self):
        J = self.num.torsional_stiffness()
        assert(SSD(self.J, J) < 0.000001)
    
    def test_shear_center(self):
        z_shear, y_shear = self.num.Shear_center()
        assert(SSD(self.ysc, y_shear)< 0.001)
        assert(SSD(self.zsc, z_shear)< 0.001)
    
    def test_all(self):
        self.test_stringer()
        self.test_crossextionalArea()
        self.test_centroid()
        self.test_momentOfinertia()
        self.test_torsional_stiffenss()
        self.test_shear_center()
        
        
run = VerifyNum().test_all()      