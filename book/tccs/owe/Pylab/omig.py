#ParametricPlot3D[{0.5 + 0.5 Sin[a Pi/180], 
#          Sin[a Pi/180], -0.5 Cos[a Pi/180], 
#	  {Thickness[0.01]}}, {a, -90, 90}, BoxRatios -> {1, 2, 1}, 
#        AxesLabel -> {"x (km)", "p (s/km)", "z (km)"}, 
#        Ticks -> {Automatic, 
#            Automatic, {0, {-0.2, "0.2"}, {-0.4, "0.4"}}}];

import matplotlib.pyplot as plt
import numpy as np
from math import pi

ax = plt.figure().add_subplot(projection='3d')

a = np.linspace(-90, 90, 200)
theta = a*pi/180
x = 0.5 + 0.5*np.sin(theta)
p = np.sin(theta)
z = -0.5*np.cos(theta)

ax.plot(x, p, z)
ax.set_xlabel('x (km)')
ax.set_ylabel('p (s/km)')
ax.set_zlabel('z (km)')
ax.set_zticks([0, -0.1, -0.2, -0.3, -0.4])
ax.set_zticklabels(["0", "0.1", "0.2", "0.3", "0.4"])

plt.savefig('junk_py.eps')
