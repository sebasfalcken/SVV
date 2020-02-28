import numpy as np
import math as m
import Properties
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


class Stress:
    def __init__(self, plane):
        self.plane = plane
        self.Izz, self.Iyy = Properties.Properties(1).MOI()
        self.z_coord, self.y_coord = Properties.Properties(1).Centroid()
        self.Ca = 0.484  # m
        self.la = 1.691  # m
        self.ha = 0.173  # m
        self.tsk = 1.1/1000  # m
        self.tsp = 2.5/1000  # m
        self.tst = 1.2/1000  # m
        self.hst = 14./1000   # m
        self.wst = 18./1000   # m
        self.nst = 13  # -
        
    def stress_zz(self, My, z):
        return My*(z-self.z_coord)/self.Iyy
    
    def stress_yy(self, Mz, y):
        return Mz*(y/self.Izz)
    
    def Print_Stress(self, My, Mz):
        #region 1
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
        stress = []
        
        for i in range(len(x)):
            s_zz = self.stress_zz(My, x[i]) 
            s_yy = self.stress_yy(Mz, y[i]) 
            stress.append(s_zz+s_yy)
        return stress, x, y
    
    def Plot(self, stress, x, y):
        
        data = np.array(stress)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        fig, axis = plt.subplots(1, 1)
        
        norm = plt.Normalize(data.min(), data.max())
        line_col = LineCollection(segments, cmap='rainbow', norm=norm)
        # Set the values used for colormapping
        line_col.set_array(data)
        line_col.set_linewidth(2)
        line = axis.add_collection(line_col)
        fig.colorbar(line, ax=axis, label=r'$\sigma_{xx} = N/m^2 $')
        
        # Use a boundary norm instead
        
        axis.set_xlim(0.05, -0.5)
        axis.set_ylim(-0.15, 0.15)
        axis.set_xlabel('z [m]')
        axis.set_title('Direct Stress distribution')
        axis.set_ylabel('y [m]')
        plt.show()

if __name__ == "__main__":
    #My = 108015.20043911
    #Mz = -29223.74268403
    
    My = 7730.1705719609545
    Mz = -37713.53236659488
    
    a = Stress(1)
    stress, x, y = a.Print_Stress(My, Mz)
    a.Plot(stress, x, y)