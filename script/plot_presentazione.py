# ::setlocal makeprg=cd\ script\ &&\ python3\ plot_presentazione.py

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

## gaussiana
x = np.linspace(-2.4, 2.4, 100)
y = np.linspace(-2.4, 2.4, 100)
x, y = np.meshgrid(x, y)
z = np.exp(-x**2 - y**2)

ax.plot_surface(x, y, z, cmap='plasma', edgecolor='black', vmin=-0.2, vmax=1.5)
    
## Path lungo
x = np.linspace(-0.04, 0.04, 100)
y = np.linspace(-2, 2, 100)
x, y = np.meshgrid(x, y)
z = np.exp(-y**2 - x**2) + 0.1
ax.plot_surface(x, y, z, color='#08d204', label='cammino più lungo')


## Path lungo dietro
y = np.linspace(0, 2, 100)
x, y = np.meshgrid(x, y)
z = np.exp(-y**2 - x**2)
#ax.plot_surface(x, y, z)

## Path corto
#for R in np.linspace(1.95, 2.05, 10):
theta = np.linspace(0, np.pi, 100)
R = 2# - 0.5 * np.sin(theta)**2
x = R*np.sin(theta)
y = R*np.cos(theta)
z = np.exp(-y**2 - x**2) + 0.1

ax.plot(x, y, z, color='red', linewidth=3, label='cammino più breve')


## Points
ax.scatter(0, 2, 0, color='orange', s=100)
ax.scatter(0, -2, 0, color='orange', s=100)


ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.view_init(elev=45, azim=-55)
ax.view_init(elev=90, azim=-90)
plt.tight_layout()
#plt.legend(loc='lower left', fontsize=15)
plt.show()

