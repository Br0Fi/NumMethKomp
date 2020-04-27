# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math
import sys, os
sys.path.append('Zettel6/')
from rk4l import *

# Konstanten fuer einheitliche Darstellung
fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Aufgabenspezifische Konstanten:
s_max = 1E6
s_min = -1E6
x0 = [0,1]
xT = 0
t0 = 0
T = 1
h = 0.001
eps = 0.001
ds = 0.001

# Hauptprogramm:

# Numerische Konstanten:




def rhs(t, x, lambwo):
    x0,x1 = x[0],x[1]
    rhs_0 = x1
    rhs_1 = - lambwo * x0
    return np.array([rhs_0,rhs_1])

def newton(x0, h, t0, T, eps, s_min, s_max, ds, rhs, sstart):
    sn = sstart  #(s_max-s_min)/2 + s_min # starting point for newton in the middle of given interval
    snp = s_max + 1000 # I'd really just need a do while
    xv0 =  x0
    while(abs(snp-sn)>eps):
        t1, x1 = rk4(t0, xv0, T, h, rhs, sn) # one shot
        xn1 = x1[-1,0] - xT # difference betweend hit and target
        t2, x2 = rk4(t0, xv0, T, h, rhs, sn+ds) # one shot
        xn2 = x2[-1,0] - xT

        deriv = (xn2 - xn1) / h
        snp = sn
        sn = sn - xn1/deriv

        if(s_min>sn or sn>s_max):
            print("Error: s_n out of bounds")
            break
    return t1,x1,sn

for sstart in [0.5, 50, 100]:
    t, x, s = newton(x0, h, t0, T, eps, s_min, s_max, ds, rhs, sstart)
    fig = plt.figure(figsize=fig_size)
    plt.plot(t,x[:,0],".")
    plt.title("$s_{start}$ = " + str(sstart) + ", $\lambda$ = " + str(round(s,3)))
    plt.grid()
    plt.ylabel('x_max')
    plt.xlabel("t_max")
    plt.tick_params(labelsize=fig_labelsize)
    plt.show()

# %%
