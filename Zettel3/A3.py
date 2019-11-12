# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math

from scipy.special import lambertw
import cmath

# Konstanten fuer einheitliche Darstellung
fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Aufgabenspezifische Konstanten:
lambwo = 2.5
y0 = 0.5

# Hauptprogramm:

# Numerische Konstanten:
t0 = -1
T = 40
h = 0.01

# Anfangswerte


# verwende "ein klassisches RK4":
def rk4_2d(t0,y0,T,h, lambwo):
    N = int((T-t0)/h)
    t = np.linspace(t0, T, N+1)
    h_inv = int(round(1/h))
    y = np.empty(N+1)
    y[:h_inv+1] = y0
    def fy(y_arg, ym1_arg):
        erg = lambwo * y_arg * (1-ym1_arg)
        return erg

    for i in range(h_inv,N):
        ky1 = fy(y[i], y[i-h_inv])
        ky2 = fy(y[i] + h/2 * ky1, y[i-h_inv])
        ky3 = fy(y[i] + h/2 * ky2, y[i-h_inv])
        ky4 = fy(y[i] + h   * ky3, y[i-h_inv])

        y[i+1] = y[i] + h/6 * (ky1 + 2*ky2 + 2*ky3 + ky4)
    return t,y

#k = 2
#N = int(T/h)
#t = np.linspace(t0, T, N+1)
#mü = lambertw(-lambwo, k)
#print(mü)
#ute = np.exp(mü*t)
#y = (1 + ute).real # not sure about how to plot the complex array...


t, y = rk4_2d(t0, y0, T, h, lambwo)

fig = plt.figure(figsize=fig_size)
plt.title("Zeitverzögerung, y0 = " +
        ", h = " + str(h) +", lambda= " + str(lambwo))
plt.plot(t, y,".")
plt.grid()
plt.ylabel('y')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
plt.savefig("Zettel3/figures/A3.png")
plt.show()
# %%

fig = plt.figure(figsize=fig_size)
plt.title("Zeitverzögerung, y0 = " +
        ", h = " + str(h) +", lambda= " + str(lambwo))
for lambwo in [1.0, 1.8, 2.8, 3.0]:
    t, y = rk4_2d(t0, y0, T, h, lambwo)
    plt.plot(t, y,".", label="lambda = " + str(lambwo))

plt.legend(prop={'size':fig_legendsize})
plt.grid()
plt.ylabel('y')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
plt.savefig("Zettel3/figures/A3c.png")
plt.show()
# %%
