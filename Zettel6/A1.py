# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math
import sys, os
sys.path.append('Zettel4/')
from rkf4 import *

# Konstanten fuer einheitliche Darstellung
fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Aufgabenspezifische Konstanten:
s_max = 0
s_min = -100
x0 = 4
xT = 1
t0 = 0
T = 1
h = 0.001
eps = 0.001
ds = 0.001

# Hauptprogramm:

# Numerische Konstanten:




def rhs(t, x):
    x0,x1 = x[0],x[1]
    rhs_0 = x1
    rhs_1 = 3/2 * x0**2
    return np.array([rhs_0,rhs_1])

def newton(x0, h, t0, T, eps, s_min, s_max, ds, rhs, sstart):
    sn = sstart  #(s_max-s_min)/2 + s_min # starting point for newton in the middle of given interval
    snp = s_max + 1000 # I'd really just need a do while
    while(abs(snp-sn)>eps):
        #nn = int(sn/h + t[0]) # position of current s_n in array
        xv0 =  [x0, sn]
        t1, x1 = rk4(t0, xv0, T, h, rhs) # one shot
        xn1 = x1[-1,0] - xT # difference betweend hit and target
        xv0 = [x0, sn+ds]
        t2, x2 = rk4(t0, xv0, T, h, rhs) # one shot
        xn2 = x2[-1,0] - xT

        deriv = (xn2 - xn1) / h
        snp = sn
        sn = sn - xn1/deriv

        if(s_min>sn or sn>s_max):
            print("Error: s_n out of bounds")
            break
    return t1,x1,sn

for sstart in [-100, -60, -20, -10, -5, 0]:
    t, x, s = newton(x0, h, t0, T, eps, s_min, s_max, ds, rhs, sstart)
    fig = plt.figure(figsize=fig_size)
    plt.plot(t,x[:,0],".")
    plt.title("$s_{start}$ = " + str(sstart) + ", $s_{opt}$ = " + str(round(s,3)))
    plt.grid()
    plt.ylabel('x_max')
    plt.xlabel("t_max")
    plt.tick_params(labelsize=fig_labelsize)
    plt.show()

# %%
