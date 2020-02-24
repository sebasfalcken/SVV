# ----------------- Imports -----------------
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
import operator

# ----------------- Mesh Grid -----------------

N_z = 81                                                     # number of chordwise stations
C_a = 0.484                                                  # [m]
i_row = np.arange(1, N_z + 2, 1)                             # ith row
theta_z = (i_row - 1) * np.pi / N_z
z = np.asarray([])                                           # intialize array for z-coordinates
for idx in range(0, N_z):
    z = np.append(z, -0.5 * (0.5 * C_a * (1 - np.cos(theta_z[idx])) + 0.5 * C_a * (1 - np.cos(theta_z[idx + 1]))))

N_x = 41                                                     # number of spanwise stations
l_a = 1.691                                                  # [m]
i_col = np.arange(1, N_x + 2, 1)                             # ith column
theta_x = (i_col - 1) * np.pi / N_x
x = np.asarray([])                                           # intialize array for x-coordinates
for idx in range(0, N_x):
    x = np.append(x, 0.5 * (0.5 * l_a * (1 - np.cos(theta_x[idx])) + 0.5 * l_a * (1 - np.cos(theta_x[idx + 1]))))

x_mesh, y_mesh = np.meshgrid(x, z)
'''
fig, ax = plt.subplots()
ax.plot(x_mesh.flatten(), y_mesh.flatten(), ".")  # ","
ax.set(xlabel='x [m]', ylabel='z [m]',
       title='Mesh Grid')
'''
# fig.clear()
# plt.close()
# plt.show()
# fig.savefig("meshgrid.png")

# ----------------- 2D Aileron Loading -----------------

file = open("aerodynamicloadcrj700.dat").read()
rows = file.split("\n")                             # spanwise stations
load_array = []

for row in rows:
    row = row.split(",")
    row = list(np.float_(row))
    load_array.append(row)


load_array = np.asarray(load_array)

load = load_array[41]
span = np.arange(0, 41, 1)
'''
fig, ax = plt.subplots()
ax.plot(x, load)
ax.set(xlabel='x [m]', ylabel='Load (kN)',
       title='2D Aileron Loading')
'''
# fig.clear()
# plt.close()
# plt.show()


# ----------------- 3D Aileron Loading -----------------


x_grid = np.zeros((81, 41))
z_grid = np.zeros((81, 41))

counter = 0
for i in x:
    x_grid[:81, counter] = i
    counter += 1

counter = 0
for i in z:
    z_grid[counter, :41] = i
    counter += 1
'''
ax = plt.axes(projection='3d')
ax.plot_surface(x_grid, z_grid, load_array, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')
ax.set_zlabel('Distributed Force $[kN/m^{2}]$')
ax.set_title('3D Aerodynamic Loading')
'''
# fig.clear()
# plt.close()
# plt.show()


# ----------------- Numerical Methods -----------------


def interpolate(x, x_values, y_values):
    def basis(i):
        p = [(x - x_values[j])/(x_values[i] - x_values[j]) for j in range(k) if j != i]
        return reduce(operator.mul, p)
    assert len(x_values) != 0 and (len(x_values) == len(y_values)), 'x and y must have the same length'
    k = len(x_values)

    return sum(basis(i) * y_values[i] for i in range(k))


def integrate(f, a, b, n):
    """Compute the Riemann sum of f(x) over an interval [a,b].

        Parameters
        ----------
        f : array
            Array including values of f(x) to be integrated
        a , b : float
            Endpoints of the interval [a,b]
        n : integer
            Number of sub-intervals of equal length in the interval [a,b]

        Returns
        -------
        float
            Approximation of the integral given by the Riemann sum.
        """

    dx = abs(b - a) / n
    f_mid = (f[:-1] + f[1:]) / 2

    return np.sum(f_mid * dx)


# ----------------- Chordwise Discretization & Integration -----------------

load_station = np.asarray([])                                       # array for resultant point force in each chord

n_chord = 1000                                                     # number of chordwise stations
z_n = np.linspace(z[0], z[-1], n_chord)                            # z values of the new mesh

for station in range(len(x)):                                       # for each spanwise station of the 41
    load_i = interpolate(z_n, z, load_array[:, station])            # load array for spanwise station
    load_station = np.append(load_station, integrate(load_i, z_n[0],  z_n[-1], n_chord))    # q(x)


# ----------------- Spanwise Discretization & Integration -----------------

n_span = 1000                                                     # number of spanwise stations
x_n = np.linspace(x[0], x[-1], n_span)                            # x values of the new mesh
load_i_n = interpolate(x_n, x, load_station)                       # q(x) for the new mesh

# for station in range(len(x)):                                      # for each spanwise station of the 41


# -----------------Testing -----------------

# z_n = np.linspace(z[0], z[-1], 1000)
# loadd = load_array[:, 20]   # 0 to 40
# load_z = interpolate(z_n, z, loadd)
#
# fig, ax = plt.subplots()
# ax.scatter(z_n, load_z)
# ax.set(xlabel='z_n [m]', ylabel='load_z (kN)',
#        title='2D Aileron Loading')

#plt.show()

# torque, find point where force acts


# ----------------- Integral Calculations -----------------

C_a = z[-1]
l_a = x[-1]

i1 = integrate(load_i_n, x[0], x[-1], n_span + 1)
#print("i1 = " + str(i1))


I1 = np.asarray([])
counter = 1
for i in x_n:
    I1 = np.append(I1, integrate(load_i_n[0:counter], 0, i, counter))
    counter += 1
i2 = integrate(I1, x[0], x[-1], n_span)
#print("i2 = " + str(i2))

I2 = np.asarray([])
counter = 1
for i in x_n:
    I2 = np.append(I2, integrate(I1[0:counter], 0, i, counter))
    counter += 1
i3 = integrate(I2, x[0], x[-1], n_span)
#print("i3 = " + str(i3))

I3 = np.asarray([])
counter = 1
for i in x_n:
    I3 = np.append(I3, integrate(I2[0:counter], 0, i, counter))
    counter += 1
i4 = integrate(I3, x[0], x[-1], n_span)
#print("i4 = " + str(i4))

I4 = np.asarray([])
counter = 1
for i in x_n:
    I4 = np.append(I4, integrate(I3[0:counter], 0, i, counter))
    counter += 1


def int1(x):

    idx = (np.abs(x_n - x)).argmin()        # index

    return I1[idx]


def int2(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return I2[idx]


def int3(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return I3[idx]


def int4(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return I4[idx]
