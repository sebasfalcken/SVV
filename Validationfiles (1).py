import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
from pylab import *
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D


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

for i in range(0, len(coord)):
    coord[i][3] = coord[i][3] - 102.5

for i in range(0, len(coord)):
    for j in range(1, 4):
        coord[i][j] = coord[i][j] / 1000

# --------------------- Delete Failed Output Lists -----------

for i in range(7, 12):
    tables.pop(i)

tables.pop()

# --------------------- Reshape output arrays ---------------
x_HL = []
defl_y_HL_bending = []
defl_y_HL_jam_bent = []
defl_y_HL_jam_straight = []
defl_z_HL_bending = []
defl_z_HL_jam_bent = []
defl_z_HL_jam_straight = []

for i in range(0, 6):
    tables[i] = np.reshape(tables[i], (-1, 6))
    tables[i] = np.delete(tables[i], 1, 1)

for i in range(0, 6, 2):
    tables[i] = np.vstack((tables[i], tables[i + 1]))

for i in range(1, 4):
    tables.pop(i)

for i in range(3, 9):
    tables[i] = np.reshape(tables[i], (-1, 5))

for i in range(0, 3):
    tables[i] = tables[i][np.argsort(tables[i][:, 0])]
    tables[i] = tables[i][:6588]

for i in range(0, len(coord)):
    if coord[i][3] == 0.0:
        x_HL.append(coord[i][1])
        defl_y_HL_bending.append(tables[3][i][3])
        defl_y_HL_jam_bent.append(tables[4][i][3])
        defl_y_HL_jam_straight.append(tables[5][i][3])
        defl_z_HL_bending.append(tables[3][i][4])
        defl_z_HL_jam_bent.append(tables[4][i][4])
        defl_z_HL_jam_straight.append(tables[5][i][4])

# --------------------- Y Delfection -------------------------
fig, ax = plt.subplots()
ax.scatter(x_HL, defl_y_HL_bending)
ax.set(xlabel='x location [m]', ylabel='Deflection in y',
       title='Validation Data Deflection Bending')

fig, ax = plt.subplots()
ax.scatter(x_HL, defl_y_HL_jam_bent)
ax.set(xlabel='x location [m]', ylabel='Deflection in y',
       title='Validation Data Deflection Jam_Bent')

fig, ax = plt.subplots()
ax.scatter(x_HL, defl_y_HL_jam_straight)
ax.set(xlabel='x location [m]', ylabel='Deflection in y',
       title='Validation Data Deflection Jam_Straight')

# --------------------- Z Delfection -------------------------

fig, ax = plt.subplots()
ax.scatter(x_HL, defl_z_HL_bending)
ax.set(xlabel='x location [m]', ylabel='Deflection in z',
       title='Validation Data Deflection Bending')

fig, ax = plt.subplots()
ax.scatter(x_HL, defl_z_HL_jam_bent)
ax.set(xlabel='x location [m]', ylabel='Deflection in z',
       title='Validation Data Deflection Jam_Bent')

fig, ax = plt.subplots()
ax.scatter(x_HL, defl_z_HL_jam_straight)
ax.set(xlabel='x location [m]', ylabel='Deflection in z',
       title='Validation Data Deflection Jam_Straight')
# plt.show()


# --------------------- 3D Heat Map -------------------------


x = np.asarray(coord[:, 1:2])
y = np.asarray(coord[:, 2:3])
z = np.asarray(coord[:, 3:4])               # z should be in m
stress = np.asarray(tables[2][:, 1:2])



# 1
# max = max(stress)
# x = x.reshape(108, 61)
# y = y.reshape(108, 61)
# z = z.reshape(108, 61)
# stress = stress.reshape(108, 61)
#
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.view_init(45,60)
# colors = cm.hsv(stress/max)
# colmap = cm.ScalarMappable(cmap=cm.hsv)
# colmap.set_array(stress)
#
# ax.plot_surface(x, y, z, facecolors=colors)
# plt.show()
#

# 2

# x = x.reshape(x.shape[0])
# y = y.reshape(y.shape[0])
# z = z.reshape(z.shape[0])
# stress = stress.reshape(stress.shape[0])
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.scatter(x, y, z, c=stress, lw=0, s=30)
# plt.show()




# 3
x = x.reshape(x.shape[0])
y = y.reshape(y.shape[0])
z = z.reshape(z.shape[0])
stress = stress.reshape(stress.shape[0])

fig = plt.figure(figsize=(8,6))

ax = fig.add_subplot(111,projection='3d')


colors = cm.hsv(stress/max(stress))

colmap = cm.ScalarMappable(cmap=cm.hsv)
colmap.set_array(stress)

yg = ax.scatter(x, y, z, c=colors, marker='o')
cb = fig.colorbar(colmap)

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m] CONVERT ME TO M')


plt.show()


