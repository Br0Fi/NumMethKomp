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
a = 0.004
b = 50.0
c = 0.75
d = 0.001
e = 100.0
f = 3.0

# Hauptprogramm:

# Numerische Konstanten:
t0 = 0
T = 400
h = 0.1
# Anfangswerte, nicht spezifiziert durch Aufgabe.
x0 = 50
y0 = 60

# verwende "ein klassisches RK4":
def rk4_2d(t0,x0,y0,T,h):
    N = int(T/h)
    t = np.linspace(t0, T, N+1)
    x = np.empty(N+1)
    y = np.empty(N+1)
    x[0] = x0
    y[0] = y0
    def fx(x_arg, y_arg):
        return a*x_arg * (b-x_arg-c*y_arg)
    def fy(x_arg, y_arg):
        return d*y_arg * (e-y_arg-f*x_arg)

    for i in range(N):
        kx1 = fx(x[i],             y[i]            )
        ky1 = fy(x[i],             y[i]            )
        kx2 = fx(x[i] + h/2 * kx1, y[i] + h/2 * ky1)
        ky2 = fy(x[i] + h/2 * kx1, y[i] + h/2 * ky1)
        kx3 = fx(x[i] + h/2 * kx2, y[i] + h/2 * ky2)
        ky3 = fy(x[i] + h/2 * kx2, y[i] + h/2 * ky2)
        kx4 = fx(x[i] + h   * kx3, y[i] + h   * ky3)
        ky4 = fy(x[i] + h   * kx3, y[i] + h   * ky3)

        x[i+1] = x[i] + h/6 * (kx1 + 2*kx2 + 2*kx3 + kx4)
        y[i+1] = y[i] + h/6 * (ky1 + 2*ky2 + 2*ky3 + ky4)
    return t,x,y


fig, ax1 = plt.subplots(figsize=fig_size)
plt.title("Wettbewerbsmodell, x0 = " + str(x0) + ", y0 = " + str(y0) + ", h = " + str(h))

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
#plt.savefig("Zettel2/figures/A1-3.png")
plt.show()


# phasenraum:
if(True):
    x1 = np.linspace(-20.0,80.0,40)
    y1 = np.linspace(-10.0,133.0,40)
    x1, y1=np.meshgrid(x1, y1)
    vx=a*x1*(b-x1-c*y1)
    vy=d*y1*(e-y1-f*x1)
    fig1 = plt.figure(figsize=fig_size)
    '''
    Plot fixpoints of the system: (0,0) ,(20,40), (0,100), (50,0)
    '''
    plt.plot(50,0,'ro',markersize=10,markeredgecolor='r',markeredgewidth=2)
    plt.plot(0,100,'ro',markersize=10,markeredgecolor='r',markeredgewidth=2)
    plt.plot(20,40,'ro', markersize=10, markeredgecolor='r',markeredgewidth=2)
    plt.plot(0,0,'ro', markersize=10, markeredgecolor='r',markeredgewidth=2)
    #plt.plot(20,40, marker='o', markersize=10, markeredgecolor='r',mfc='None',markeredgewidth=2)
    #plt.plot(0,0,marker='o',markersize=10, markeredgecolor='r',mfc='None',markeredgewidth=2)
    plt.xlabel("x",fontsize=16)
    plt.ylabel("y",fontsize=16)
    plt.tick_params(labelsize=fig_labelsize)
    plt.xlim([-20,60])
    plt.ylim([-10,133])
    plt.streamplot(x1,y1,vx,vy,color='b',linewidth=2)

    #plt.quiver(x1, y1, vx, vy, pivot='middle', headwidth=3, headlength=6)

    plt.show()


# %%
