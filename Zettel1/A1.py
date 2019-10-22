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

import periodic


# Konstanten fuer einheitliche Darstellung

fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
#colors

# weitere Werte, Konstanten
# Werte von https://physics.nist.gov/cuu/Constants/index.html[0]

c = 299792458 # m/s
Ryd = 13.605693122 # eV [https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev]
K = 273.15 # kelvin
g = 9.81 # m/s^2
rad = 360 / 2 / math.pi
grad = 1/rad


def ff(t, x):
    return t-x**2
    
def euler(aa, bb, hh, x_0, ff, tt, meshpoints):
    
    xx = np.zeros(meshpoints)
    xx[0] = x_0
    for i in range(1, meshpoints):
        xx[i] = xx[i-1] + hh * ff(tt[i], xx[i-1])
    return xx
    
aa = 0.0
bb = 9.0
hh = 0.05
x_0 = [-0.7, 0.0, 1.0, 3.0]

fig = plt.figure(figsize=fig_size)
for aufg in range(4):
    meshpoints = int(round((bb-aa)/hh+0.5)) # aufrunden
    tt = np.arange(meshpoints)
    tt = aa + tt * hh
    xx = euler(aa, bb, hh, x_0[aufg], ff,tt,meshpoints)
    plt.plot(tt, xx,".",label="x_0= " + str(x_0[aufg]))

plt.grid()
plt.legend(prop={'size':fig_legendsize})
plt.ylabel('x')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
plt.savefig("Zettel1/figures/A1-a.png")
plt.show()
# %%

aa = 0.0
bb = 900.0
hh = [0.05, 0.025]
x_0 = 0.0

fig = plt.figure(figsize=fig_size)
for aufg in range(2):
    meshpoints = int(round((bb-aa)/hh[aufg]+0.5)) # aufrunden
    tt = np.arange(meshpoints)
    tt = aa + tt * hh[aufg]
    xx = euler(aa, bb, hh[aufg], x_0, ff,tt,meshpoints)
    plt.plot(tt, xx,".",label="h= " + str(hh[aufg]))

plt.grid()
plt.legend(prop={'size':fig_legendsize})
plt.ylabel('x')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
plt.savefig("Zettel1/figures/A1-bc.png")
plt.show()
print('Interpretation: Halbierung der Schrittweite kann den Unterschied\n'
    'zwischen Chaos und Lösung darstellen. Bei zu kleiner Schrittweite\n' 
    'findet eine Oszillation zwischen Funktionswerten statt.')
# %%
