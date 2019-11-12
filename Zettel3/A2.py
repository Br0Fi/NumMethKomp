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
sigma = 10
bb = 8/3
rr = 0.5
part = "a"

# Hauptprogramm:

# Numerische Konstanten:
t0 = 0
T = 200
h = 0.005
# Anfangswerte
x0 = 50
y0 = 50
z0 = 50
low_cut = 0
xy_cutoff = 0 + low_cut

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
    def fz(x_arg, y_arg, t_arg):
        return x_arg * y_arg - bb * z_arg

    for i in range(N):
        kx1 = fx(x[i],             y[i]            , z[i], t[i] +h/2)
        ky1 = fy(x[i],             y[i]            , z[i], t[i] +h/2)
        kz1 = fz(x[i],             y[i]            , z[i], t[i] +h/2)
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


fig, ax1 = plt.subplots(figsize=fig_size)
plt.title("Lorenz-System, x0 = " + str(round(x0, 5)) +", y0 = " + str(y0) +
        ", z0 = " + str(z0) + ", h = " + str(h) +", r= " + str(rr))

t, x, y, z = rk4_2d(t0, x0, y0, z0, T, h)

color = "tab:red"
plt.legend(prop={'size':fig_legendsize})
ax1.set_xlabel('t')
ax1.set_ylabel('x, z', color=color)
ax1.plot(t[low_cut:xy_cutoff], x[low_cut:xy_cutoff], label="x(t)", color=color)
ax1.plot(t[low_cut:xy_cutoff], x[low_cut:xy_cutoff], label="z(t)", color="tab:green")
plt.grid()
ax1.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis
color = "tab:blue"
ax2.set_ylabel('y', color=color)
ax2.plot(t[low_cut:xy_cutoff],y[low_cut:xy_cutoff], label="y(t)", color=color)
ax2.tick_params(axis='y', labelsize=fig_labelsize, labelcolor=color)
fig.tight_layout()
#plt.savefig("Zettel3/figures/A2-" + part + ".png")
plt.show()

# %%

# Return Map:
def return_map(t,x):
    N = int(T/h)
    xmax = x[(x > np.roll(x,1)) & (x > np.roll(x,-1))]
    tmax = t[(x > np.roll(x,1)) & (x > np.roll(x,-1))]

    return tmax, xmax

fig = plt.figure(figsize=fig_size)
tmax, zmax = return_map(t,x)
plt.plot(tmax, zmax,".")
plt.grid()
plt.ylabel('z_max')
plt.xlabel("t_max")
plt.tick_params(labelsize=fig_labelsize)
plt.show()

fig = plt.figure(figsize=fig_size)
plt.title("Return Plot")
#plt.xlim([2,2.4])
#plt.ylim([2,2.4])
zmax = zmax[int(len(zmax)/2):-1] # crop tune in time and crop last maximum (not actual maximum)
phi_n = zmax[:-1]
phi_n1 = zmax[1:]
plt.plot(phi_n, phi_n1,".")

plt.grid()
plt.ylabel('z_max(n+1)')
plt.xlabel("Z_max(n)")
plt.tick_params(labelsize=fig_labelsize)
#plt.savefig("Zettel3/figures/A2-" + part + "-return.png")
plt.show()
# %%
