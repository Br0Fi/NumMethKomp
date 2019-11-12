# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math

# Konstanten fuer einheitliche Darstellung
fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Aufgabenspezifische Konstanten:
a_0 = 2
eps = 0.3
b = 0.02
c = 0.0002
d = 0.8

# Hauptprogramm:

# Numerische Konstanten:
t0 = 0
T = 50
h = 0.01
# Anfangswerte
x0 = 5000
y0 = 120

# verwende "ein klassisches RK4":
def rk4_2d(t0,x0,y0,T,h):
    N = int(T/h)
    t = np.linspace(t0, T, N+1)
    x = np.empty(N+1)
    y = np.empty(N+1)
    x[0] = x0
    y[0] = y0
    def fx(x_arg, y_arg, t_arg):
        a = a_0* (1 + eps * math.sin(math.pi * t_arg))
        return a*x_arg - b*x_arg*y_arg
    def fy(x_arg, y_arg, t_arg):
        return c*x_arg*y_arg - d*y_arg

    for i in range(N):
        kx1 = fx(x[i],             y[i]            , t[i] +h/2)
        ky1 = fy(x[i],             y[i]            , t[i] +h/2)
        kx2 = fx(x[i] + h/2 * kx1, y[i] + h/2 * ky1, t[i] +h/2)
        ky2 = fy(x[i] + h/2 * kx1, y[i] + h/2 * ky1, t[i] +h/2)
        kx3 = fx(x[i] + h/2 * kx2, y[i] + h/2 * ky2, t[i] +h/2)
        ky3 = fy(x[i] + h/2 * kx2, y[i] + h/2 * ky2, t[i] +h/2)
        kx4 = fx(x[i] + h   * kx3, y[i] + h   * ky3, t[i] +h)
        ky4 = fy(x[i] + h   * kx3, y[i] + h   * ky3, t[i] +h)

        x[i+1] = x[i] + h/6 * (kx1 + 2*kx2 + 2*kx3 + kx4)
        y[i+1] = y[i] + h/6 * (ky1 + 2*ky2 + 2*ky3 + ky4)
    return t,x,y


fig, ax1 = plt.subplots(figsize=fig_size)
plt.title("RB-Modell, x0 = " + str(x0) + ", y0 = " + str(y0) + ", h = " + str(h) +", epsilon= " + str(eps))

t, x, y = rk4_2d(t0, x0, y0, T, h)

color = "tab:red"
ax1.set_xlabel('t')
ax1.set_ylabel('x', color=color)
ax1.plot(t, x, label="x(t)", color=color)
plt.grid()
plt.legend(prop={'size':fig_legendsize})
ax1.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis
color = "tab:blue"
ax2.set_ylabel('y', color=color)
ax2.plot(t,y, label="y(t)", color=color) # not sure how to display this label.
ax2.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
fig.tight_layout()
#plt.savefig("Zettel2/figures/A2-2-2.png")
plt.show()


# %%
