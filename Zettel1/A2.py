# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import statistics as stat
import scipy as sci
import scipy.fftpack
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math
from scipy import optimize
import glob
import os


# Konstanten fuer einheitliche Darstellung

fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

def euler(aa, bb, hh, x_0, ff, tt, meshpoints):

    xx = np.zeros(meshpoints)
    xx[0] = x_0
    for i in range(1, meshpoints):
        xx[i] = xx[i-1] + hh * ff(tt[i-1], xx[i-1])
    return xx
def heun(aa, bb, hh, x_0, ff, tt, meshpoints):

    xx = np.zeros(meshpoints)
    xx[0] = x_0
    for i in range(1, meshpoints):
        xx[i] = xx[i-1] + hh * ff(tt[i-1], xx[i-1])
        xx[i] = xx[i-1] + hh/2 * (ff(tt[i-1], xx[i-1]) + ff(tt[i], xx[i]))
    return xx

v_0 = 1
omega = 1
x_0 = 0
y_0 =v_0

x = np.linspace (-2*np.pi,2*np.pi,30)
y = np.linspace (-2*np.pi,2*np.pi,30)
x_grid,y_grid = np.meshgrid(x,y)
vx = y_grid
vy= -omega**2*x_grid
Z = 0.5*y_grid**2 + 0.5 * omega**2 * x_grid**2

fig, ax = plt.subplots(figsize=fig_size)
plt.quiver(x_grid,y_grid,vx,vy)
CS = plt.contour(x_grid, y_grid, Z)
ax.clabel(CS, inline=1, fontsize=10)
plt.ylabel('y')
plt.xlabel("x")
ax.set_title('Phasenraum')
plt.savefig("Zettel1/figures/A2-phasenraum.png")
plt.show()
# %%

aa = 0
bb = 20 * math.pi
hh = [0.05, 0.025, 0.001]
def ff(t, x):
    return -1 * omega**2 * x * t + v_0

for aufg in range(3):
    meshpoints = int(round((bb-aa)/hh[aufg]+0.5)) # aufrunden
    tt = np.arange(meshpoints)
    tt = aa + tt * hh[aufg]
    xx = euler(aa, bb, hh[aufg], x_0, ff,tt,meshpoints)
    fig = plt.figure(figsize=fig_size)
    plt.ylim(bottom = -0.1, top = 1.0)
    plt.plot(tt, xx,".",label="h= " + str(hh[aufg]))
    plt.grid()
    plt.legend(prop={'size':fig_legendsize})
    plt.ylabel('x')
    plt.xlabel("t")
    plt.tick_params(labelsize=fig_labelsize)
    fig.suptitle('Euler-Verfahren')
    plt.savefig("Zettel1/figures/A2-euler-" + str(aufg+1) + ".png")
    plt.show()
# %%
hh = 0.05
meshpoints = int(round((bb-aa)/hh+0.5)) # aufrunden
tt = np.arange(meshpoints)
tt = aa + tt * hh
xx = heun(aa, bb, hh, x_0, ff,tt,meshpoints)
fig = plt.figure(figsize=fig_size)
plt.ylim(bottom = -0.1, top = 1.2)
plt.plot(tt, xx,".",label="h= " + str(hh))
plt.grid()
plt.legend(prop={'size':fig_legendsize})
plt.ylabel('x')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
fig.suptitle('Heun-Verfahren')
plt.savefig("Zettel1/figures/A2-heun.png")
plt.show()
# %%
