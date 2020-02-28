from Properties import Properties
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
import operator
from mpl_toolkits import mplot3d
import plotly.graph_objects as go

# ------------------------- Open B737.rpt -----------------------

f = open("B737.rpt").read()

# ------------------------- Separate tables ----------------------

tables = f.split("-----------------------\n")

# ---------------------- Cull Headers & Footers ------------------

for i in range(0, len(tables) - 1, 2):
    tables[i] = tables[i].split()
    tables[i] = tables[i][:-88]

for i in range(1, len(tables) - 1, 2):
    tables[i] = tables[i].split()
    tables[i] = tables[i][:-50]

tables[-1] = tables[-1].split()
tables[-1] = tables[-1][:-27]

tables = tables[1:19]

# ----------------------Setup arrays----------------------------

for i in range(0, 6):
    tables[i] = tables[i][:tables[i].index('Minimum')]

for i in range(0, len(tables)):
    for j in range(0, len(tables[i])):
        tables[i][j] = float(tables[i][j])

# -------------------- Open B737.inp --------------------------

f1 = open("B737.inp").read()

# -------------------- Cull Header and Footer -----------------

coord = f1.split("*Part, name=B737\n*Node")
coord = coord[1]

coord = coord.split()
coord = coord[:26352]

# -------------------- Setup Coordinate Array -----------------

for i in range(0, len(coord)):
    coord[i] = coord[i].replace(",", "")
    coord[i] = float(coord[i])

coord = np.reshape(coord, (-1, 4))

# --------------------- Delete Failed Output Lists -----------

for i in range(7, 12):
    tables.pop(i)

tables.pop()

tables[0] = np.reshape(tables[0], (-1, 6))
tables[0] = np.delete(tables[0], 1, 1)

# --------------------- Reshape Output Arrays -----------
for i in range(0,6):
    tables[i] = np.reshape(tables[i],(-1,6))
    tables[i] = np.delete(tables[i],1,1)


for i in range(0,6,2):
    tables[i] = np.vstack((tables[i],tables[i+1]))

for i in range(1,4):
    tables.pop(i)

for i in range(3,9):
    tables[i] = np.reshape(tables[i],(-1,5))

for i in range(0,3):
    tables[i] = tables[i][np.argsort(tables[i][:,0])]
    tables[i] = tables[i][:6588]

# --------------------- Convert to Dictionary -----------

data = {}
for i in range(len(coord)):
    data.update({coord[i][0]: coord[i][1:4]})

# for element in data:
#     element.append()


""" 
9 Tables 
Table 0             Bending             Stresses                Region 1 & 2               6588
Table 1             Jam_Bent            Stresses                Region 1  & 2              6588
Table 2             Jam_Straight        Stresses                Region 1 & 2               6588
Table 3             
Table 6             Bending             Deflections             Region 1                6634




"""
x = np.asarray(coord[: , 1:2])
y = np.asarray(coord[: , 2:3])
z = np.asarray(coord[: , 3:4])



lst = np.array(([21, 22, 23], [11, 22, 33], [43, 77, 89]))
# print(lst[: , 0:1])


# ax = plt.axes(projection='3d')
# ax.plot_surface(x, y, z, rstride=1, cstride=1,
#                 cmap='viridis', edgecolor='none');
# plt.show()
#


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# ax.scatter(x, y, z, c='y', marker='o')
#
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
#
# plt.show()


xx, yy= np.meshgrid(np.unique(x), np.unique(y))
x_vals, x_idx = np.unique(x, return_inverse=True)
y_vals, y_idx = np.unique(y, return_inverse=True)
vals_array = np.empty(x_vals.shape + y_vals.shape)
print(len(z))
print(len(x_vals))
print(len(y_vals))


vals_array.fill(np.nan) # or whatever your desired missing data flag is
print(np.shape(vals_array))

vals_array[x_idx, y_idx] = z
zz = vals_array.T

