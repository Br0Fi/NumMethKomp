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
ks = 0.799
ds = 9.44
ns = 14.68
fs = 2.1
omega = 2.5

# Hauptprogramm:

# Numerische Konstanten:
t0 = 0
T = 50
h = 0.001
# Anfangswerte
x0 = -1.5*math.pi/2
y0 = 0

# verwende "ein klassisches RK4":
def rk4_2d(t0,x0,y0,T,h):
    N = int(T/h)
    t = np.linspace(t0, T, N+1)
    x = np.empty(N+1)
    y = np.empty(N+1)
    x[0] = x0
    y[0] = y0
    def fx(x_arg, y_arg, t_arg):
        return y_arg
    def fy(x_arg, y_arg, t_arg):
        return -ds*x_arg + fs*math.sin(omega*t_arg) \
               - ks*y_arg + ns*math.sin(x_arg)

    for i in range(N):
        kx1 = fx(x[i],             y[i]            , t[i] +h/2)
        ky1 = fy(x[i],             y[i]            , t[i] +h/2)
        kx2 = fx(x[i] + h/2 * kx1, y[i] + h/2 * ky1, t[i] +h/2)
        ky2 = fy(x[i] + h/2 * kx1, y[i] + h/2 * ky1, t[i] +h/2)
        kx3 = fx(x[i] + h/2 * kx2, y[i] + h/2 * ky2, t[i] +h/2)
        ky3 = fy(x[i] + h/2 * kx2, y[i] + h/2 * ky2, t[i] +h/2)
        kx4 = fx(x[i] + h   * kx3, y[i] + h   * ky3, t[i] +h/2)
        ky4 = fy(x[i] + h   * kx3, y[i] + h   * ky3, t[i] +h/2)

        x[i+1] = x[i] + h/6 * (kx1 + 2*kx2 + 2*kx3 + kx4)
        y[i+1] = y[i] + h/6 * (ky1 + 2*ky2 + 2*ky3 + ky4)
    return t,x,y


fig, ax1 = plt.subplots(figsize=fig_size)
plt.title("Pohlsches Rad, x0 = " + str(round(x0, 5)) +
        ", y0 = " + str(y0) + ", h = " + str(h) +", Omega= " + str(omega))

t, x, y = rk4_2d(t0, x0, y0, T, h)

color = "tab:red"
ax1.set_xlabel('t')
ax1.set_ylabel('x (grad)', color=color)
ax1.plot(t, x, label="x(t)", color=color)
plt.grid()
#plt.legend(prop={'size':fig_legendsize})
ax1.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis
color = "tab:blue"
ax2.set_ylabel('y', color=color)
ax2.plot(t,y, label="y(t)", color=color) # not sure how to display this label.
ax2.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
fig.tight_layout()
#plt.savefig("Zettel3/figures/A1-a.png")
plt.show()

# %%

# Return Map:
def return_map(t,x):
    N = int(T/h)
    xmax = x[(x > np.roll(x,1)) & (x > np.roll(x,-1))]
    tmax = t[(x > np.roll(x,1)) & (x > np.roll(x,-1))]
    
    return tmax, xmax

fig = plt.figure(figsize=fig_size)
tmax, xmax = return_map(t,x)
plt.plot(tmax, xmax,".")
plt.grid()
plt.ylabel('x_max')
plt.xlabel("t_max")
plt.tick_params(labelsize=fig_labelsize)
plt.show()

fig = plt.figure(figsize=fig_size)
plt.title("Return Plot")
plt.xlim([-3,3])
plt.ylim([-3,3])
xmax = xmax[int(len(xmax)/2):] # crop tune in time
phi_n = xmax[:-1]
phi_n1 = xmax[1:]
plt.plot(phi_n, phi_n1,".")

plt.grid()
plt.ylabel('phi_(n+1)')
plt.xlabel("phi_n")
plt.tick_params(labelsize=fig_labelsize)
#plt.savefig("Zettel3/figures/A1-a_return.png")
plt.show()
# %%
