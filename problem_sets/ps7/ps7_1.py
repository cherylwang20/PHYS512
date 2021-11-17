import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

dat = np.loadtxt('rand_points.txt')
dat_py = np.loadtxt('rand_points_py.txt')
x = dat[:,0]; x_py = dat_py[:,0];
y = dat[:,1]; y_py = dat_py[:,1];
z = dat[:,2]; z_py = dat_py[:,2]



fig = plt.figure(figsize = (8, 8))
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, 'gray', s = 0.05)
plt.show()

fig = plt.figure(figsize = (8, 8))
ax = plt.axes(projection='3d')
ax.scatter3D(x_py, y_py, z_py, 'gray', s = 0.05)
plt.show()