# Please do not change anything in code before getting back to me
# ----------------- Imports -----------------
from Properties import Properties
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
import operator
from mpl_toolkits import mplot3d

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

fig, ax = plt.subplots()
ax.plot(x_mesh.flatten(), y_mesh.flatten(), ".")  # ","
ax.set(xlabel='x [m]', ylabel='z [m]',
       title='Mesh Grid')

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

fig, ax = plt.subplots()
ax.plot(x, load)
ax.set(xlabel='x [m]', ylabel='Load (kN)',
       title='2D Aileron Loading')

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

ax = plt.axes(projection='3d')
ax.plot_surface(x_grid, z_grid, load_array, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')
ax.set_zlabel('Distributed Force $[kN/m^{2}]$')
ax.set_title('3D Aerodynamic Loading')

# ----------------- Numerical Methods -----------------


def interpolate(x, x_values, y_values):
    """Interpolates a function y(x) given a discrete dataset x_values and y_values.

            Parameters
            ----------
            x : array or float
                Interpolation mesh. Array including all values you would like y(x) to be evaluated at

            x_values : numpy array
                Array including all x values in the dataset.

            y_values : numpy array
                Array including all y values in the dataset.

            Returns
            -------
            numpy array
                Array including values of y(x) for each x_i in the array x
            """

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

n_chord = 1000                                                      # number of chordwise stations
z_n = np.linspace(z[0], z[-1], n_chord)                             # z values of the new mesh

for station in range(len(x)):                                       # for each spanwise station of the 41
    load_i = interpolate(z_n, z, load_array[:, station])            # load array for spanwise station
    load_station = np.append(load_station, integrate(load_i, z_n[0],  z_n[-1], n_chord))    # q(x) tilde


# ----------------- Spanwise Discretization -----------------

n_span = 1000                                                       # number of spanwise stations
x_n = np.linspace(x[0], x[-1], n_span)                              # x values of the new mesh
load_i_n = interpolate(x_n, x, load_station)                        # q(x) for the new mesh

# ----------------- Spanwise Integral Calculations -----------------

C_a = z[-1]
l_a = x[-1]
z_sc, x_sc = Properties(1).Shear_center()


Q1 = np.asarray([])
counter = 1
for i in x_n:
    Q1 = np.append(Q1, integrate(load_i_n[0:counter], 0, i, counter))
    counter += 1

Q2 = np.asarray([])
counter = 1
for i in x_n:
    Q2 = np.append(Q2, integrate(Q1[0:counter], 0, i, counter))
    counter += 1

Q3 = np.asarray([])
counter = 1
for i in x_n:
    Q3 = np.append(Q3, integrate(Q2[0:counter], 0, i, counter))
    counter += 1

Q4 = np.asarray([])
counter = 1
for i in x_n:
    Q4 = np.append(Q4, integrate(Q3[0:counter], 0, i, counter))
    counter += 1


Tau = np.asarray([])
for station in range(len(x)):                                       # for each spanwise station of the 41
    load_i = interpolate(z_n, z, load_array[:, station])            # load array for spanwise station
    Tau = np.append(Tau, integrate(load_i * (z_n - z_sc), z_n[0],  z_n[-1], n_chord))    # q(x) tilde * (z - z tilde)


Tau1 = np.asarray([])
counter = 1
for i in x_n:
    Tau1 = np.append(Tau1, integrate(Tau[0:counter], 0, i, counter))
    counter += 1

Tau2 = np.asarray([])
for i in x_n:
    Tau2 = np.append(Tau2, integrate(Tau1[0:counter], 0, i, counter))
    counter += 1


def q1(x):

    idx = (np.abs(x_n - x)).argmin()        # index

    return Q1[idx]


def q2(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return Q2[idx]


def q3(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return Q3[idx]


def q4(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return Q4[idx]


def tau1(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return Tau1[idx]


def tau2(x):
    idx = (np.abs(x_n - x)).argmin()        # index

    return Tau2[idx]


# -----------------Testing -----------------
plt.close('all')
# z_n = np.linspace(z[0], z[-1], 1000)
# loadd = load_array[:, 20]   # 0 to 40
# load_z = interpolate(z_n, z, loadd)
#
# fig, ax = plt.subplots()
# ax.scatter(z_n, load_z)
# ax.set(xlabel='z_n [m]', ylabel='load_z (kN)',
#        title='2D Aileron Loading')

# plt.show()                                                   PLOTS ALL FIGURES

# Please do not change anything in code before getting back to me
