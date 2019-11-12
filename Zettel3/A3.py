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
T = 25
h = 0.002
# Anfangswerte


# verwende "ein klassisches RK4":
def rk4_2d(t0,x0,y0,z0,T,h):
    N = int(T/h)
    t = np.linspace(t0, T, N+1)
    x = np.empty(N+1)
    y = np.empty(N+1)
    z = np.empty(N+1)
    x[0] = x0
    y[0] = y0
    z[0] = y0
    def fx(x_arg, y_arg, z_arg, t_arg):
        return sigma * (y_arg - x_arg)
    def fy(x_arg, y_arg, z_arg, t_arg):
        return rr*x_arg - y_arg - x_arg * z_arg
    def fz(x_arg, y_arg, z_arg, t_arg):
        return x_arg * y_arg - bb * z_arg

    for i in range(N):
        kx1 = fx(x[i],             y[i]            , z[i], t[i])
        ky1 = fy(x[i],             y[i]            , z[i], t[i])
        kz1 = fz(x[i],             y[i]            , z[i], t[i])
        kx2 = fx(x[i] + h/2 * kx1, y[i] + h/2 * ky1, z[i]+h/2*kz1, t[i] +h/2)
        ky2 = fy(x[i] + h/2 * kx1, y[i] + h/2 * ky1, z[i]+h/2*kz1, t[i] +h/2)
        kz2 = fz(x[i] + h/2 * kx1, y[i] + h/2 * ky1, z[i]+h/2*kz1, t[i] +h/2)
        kx3 = fx(x[i] + h/2 * kx2, y[i] + h/2 * ky2, z[i]+h/2*kz2, t[i] +h/2)
        ky3 = fy(x[i] + h/2 * kx2, y[i] + h/2 * ky2, z[i]+h/2*kz2, t[i] +h/2)
        kz3 = fz(x[i] + h/2 * kx2, y[i] + h/2 * ky2, z[i]+h/2*kz2, t[i] +h/2)
        kx4 = fx(x[i] + h   * kx3, y[i] + h   * ky3, z[i]+h*kz3, t[i] +h)
        ky4 = fy(x[i] + h   * kx3, y[i] + h   * ky3, z[i]+h*kz3, t[i] +h)
        kz4 = fz(x[i] + h   * kx3, y[i] + h   * ky3, z[i]+h*kz3, t[i] +h)

        x[i+1] = x[i] + h/6 * (kx1 + 2*kx2 + 2*kx3 + kx4)
        y[i+1] = y[i] + h/6 * (ky1 + 2*ky2 + 2*ky3 + ky4)
        z[i+1] = z[i] + h/6 * (kz1 + 2*kz2 + 2*kz3 + kz4)
    return t,x,y,z

k = 2
N = int(T/h)
t = np.linspace(t0, T, N+1)
mü = lambertw(-lambwo, k)
print(mü)
ute = np.exp(mü*t)
y = (1 + ute).real # not sure about how to plot the complex array...

fig = plt.figure(figsize=fig_size)
plt.plot(t, y,".",label="k= " + str(k))
plt.grid()
plt.legend(prop={'size':fig_legendsize})
plt.ylabel('y')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
plt.title('Zeitverzögerung')
#plt.savefig("Zettel3/figures/A3.png")
plt.show()
# %%
