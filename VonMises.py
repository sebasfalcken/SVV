import numpy as np
import math as m
import Direct_stress
import Shearflow
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

class Von_Mises:
    def __init__(self, My, Mz, Sy, Sz, T):
        self.My = My
        self.Mz = Mz
        self.Sy = Sy
        self.Sz = Sz
        self.T = T
        self.tsk = 1.1/1000  # m
        self.tsp = 2.5/1000
    
    def Print_VonMises(self, data):
        stress, x, y = Direct_stress.Stress(1).Print_Stress(self.My, self.Mz)
        #data, x, y = Shearflow.ShearFlow(1).Print_shear(self.Sy, self.Sz, self.T)
        VonMises = []
        max_shear = 0
        for i in range(len(stress)):
            if i>= 99 and i<= 199:
                th = self.tsp
            if i>= 399 and i<= 499:
                th = self.tsp
            else:
                th = self.tsk
            if (abs(data[i]/th)) > max_shear:
                max_shear = abs(data[i]/th)
            VonMises.append(m.sqrt((stress[i]**2) + 3*(data[i]/th)**2))  
        print(max_shear)
        return np.array(VonMises), x, y
    
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
        fig.colorbar(line, ax=axis, label=r'$\sigma_{xx} = N/m^2 $')
        
        # Use a boundary norm instead
        
        axis.set_xlim(0.05, -0.5)
        axis.set_ylim(-0.15, 0.15)
        axis.set_xlabel('z [m]')
        axis.set_title('Von Mises stress distribution')
        axis.set_ylabel('y [m]')
        plt.show()    
        
if __name__ == "__main__":
    #My = 108015.20043911
    #Mz = -29223.74268403
    #Sy = -4985.0292997
    #Sz = -73359.55813566
    #T = -1743.7425962
    
    My = 51898.94614653427
    Mz = -11162.435493672943
    T = -99.59640250928304
    Sy = 25930.735282179103
    Sz = -117684.68513953342
    
    b = Direct_stress.Stress(1)
    stress, x, y = b.Print_Stress(My,Mz)
    b.Plot(stress, x, y)
    
    c = Shearflow.ShearFlow(1)
    data, x, y = c.Print_shear(Sy, Sz, T)
    c.Plot(data, x, y)
    
    a = Von_Mises(My, Mz, Sy, Sz, T)
    data, x, y = a.Print_VonMises(data)
    a.Plot(data, x, y)