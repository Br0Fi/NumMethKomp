# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math

# Spezifisch für die Art, wie ich python ausführe. Funktioniert aber vmtl. auch für Ausführen von der Kommandozeile.
import sys, os
sys.path.append('Zettel4/')
from rkf4 import *

# Parameter fuer einheitliche Darstellung
fig_size = (10, 6)
fig_fontsize = 18
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Numerische Konstanten:
t0 = 0.
N = 2**13
h_start = 0.1
eps_0 = 1E-4
beta_ctrl = 0.9
# Anfangswerte
x0 = np.array([2.0,0.0])

drop = 1/4. # Anteil der Lösung, die nicht für die Plots berücksichtigt werden soll, damit Probleme mit dem transienten Zustand (Einschwingvorgang) am Anfang der Simulation vermieden werden

## RIGHT-HAND SIDE OF ODE:
def rhs(t, x):
    x0,x1 = x[0],x[1]
    rhs0 = x1
    rhs1 = -mu*(x0**2-1)*x1-x0
    return np.array([rhs0, rhs1])

mu_values = np.array([0.5, 3.0, 10.0, 20.0])
teil = ['a)', 'b)', 'c)', 'd)']

for i in range(4):
    mu = mu_values[i]

    t, x = rkf(t0, x0, N, h_start, rhs, eps_0, beta_ctrl)

    #vergiss den Einschwingvorgang am Anfang
    t = t[int(len(t)*drop):]
    x = x[int(len(x)*drop):]

    #plotte Phasenraum
    plt.figure(figsize=fig_size)
    plt.title("Van-der-Pol-Oszi, $\mu$ = {:.3g}".format(mu), fontsize=fig_fontsize)
    plt.xlabel('$x$',fontsize=fig_fontsize)
    plt.ylabel('$y$',fontsize=fig_fontsize)

    plt.plot(x[:,0],x[:,1],'g',linewidth=2)

    #plt.savefig('B3_A1{}.png'.format(teil[i]))

plt.show()
