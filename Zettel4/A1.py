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
from rk4 import *

# Parameter fuer einheitliche Darstellung
fig_size = (10, 6)
fig_fontsize = 18
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Numerische Konstanten:
t0 = 0.
T = 400.
h = 0.01
# Anfangswerte
x0 = np.array([2.0,0.0])

drop = 1/4. # Anteil der Lösung, die nicht für die Plots berücksichtigt werden soll, damit Probleme mit dem transienten Zustand (Einschwingvorgang) am Anfang der Simulation vermieden werden

## RIGHT-HAND SIDE OF ODE:
def rhs(t, x):
    x0,x1 = x[0],x[1]
    rhs0 = x1
    rhs1 = -mu*(x0**2-1)*x1-x0
    return np.array([rhs0, rhs1])

## APPLY RK4-METHOD AND DO PLOTS FOR DIFFERENT VALUES OF mu:

mu_values = np.array([0.5, 3.0, 10.0, 20.0])
teil = ['a)', 'b)', 'c)', 'd)']

for i in range(4):
    mu = mu_values[i]

    t, x = rk4(t0, x0, T, h, rhs)

    #vergiss den Einschwingvorgang am Anfang
    t = t[int(len(t)*drop):]
    x = x[int(len(x)*drop):]

    #plotte Phasenraum
    plt.figure(figsize=fig_size)
    plt.title("Pohlsches Rad, $\mu$ = {:.3g}".format(mu), fontsize=fig_fontsize)
    plt.xlabel('$x$',fontsize=fig_fontsize)
    plt.ylabel('$y$',fontsize=fig_fontsize)

    plt.plot(x[:,0],x[:,1],'g',linewidth=2)

    #plt.savefig('B3_A1{}.png'.format(teil[i]))

plt.show()
