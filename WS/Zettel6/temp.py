# Numerische Methoden fuer komplexe Systeme I
# Blatt 2, Aufgabe 1
# Lösung des AWP x'=ax(b-x-cy), y'=dy(e-y-fx), x(t_0)=x0, y(0)=y0 over [t_0,T] via RK4-Methode mit Schrittweite h

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('Zettel4/')
from rkf4 import *

## KONSTANTEN:
a = 0.004
b = 50.
c = 0.75
d = 0.001
e = 100.
f = 3.

## RIGHT-HAND SIDE OF ODE:
def rhs(t, x):
    x0,x1 = x[0],x[1]
    rhs_0 = a*x0*(b-x0-c*x1)
    rhs_1 = d*x1*(e-x1-f*x0)
    return np.array([rhs_0,rhs_1])

## PLOT:
def plot(t0, x0, T, h, rhs):
    t, x = rk4(t0, x0, T, h, rhs)
    plt.subplot(121)
    plt.plot(x[:,0],x[:,1],linewidth=3)
    plt.subplot(122)
    plt.plot(t,x[:,0],linewidth=2)
    plt.plot(t,x[:,1],linewidth=2)

## APPLY RK4-METHOD AND DO PLOTS FOR DIFFERENT INITIAL CONDITIONS:
t0 = 0.
T = 200.
h = 0.01

plt.figure('Blatt 2, Aufgabe 1')

plt.subplot(121)
plt.xlabel("x",fontsize=16)
plt.ylabel("y",fontsize=16)
plt.xlim([-5,80])
plt.ylim([-5,105])
plt.grid()
plt.plot(50,0,marker='o',markersize=10,markeredgecolor='r',markeredgewidth=2,markerfacecolor='r')
plt.plot(0,100,marker='o',markersize=10,markeredgecolor='r',markeredgewidth=2,markerfacecolor='r')
plt.plot(20,40, marker='o', markersize=10, markeredgecolor='r',markeredgewidth=2,markerfacecolor='None')
plt.plot(0,0,marker='o',markersize=10, markeredgecolor='r',markeredgewidth=2,markerfacecolor='None')

plt.subplot(122)
plt.xlabel("t",fontsize=16)
plt.ylabel("x,y",fontsize=16)
plt.xlim([t0,T])
plt.ylim([-1,105])
plt.grid()

plot(t0,np.array([70.,100.]),T,h, rhs)
plot(t0,np.array([70.,20.]),T,h, rhs)
plot(t0,np.array([10.,10.]),T,h, rhs)
plot(t0,np.array([50.,70.]),T,h, rhs)
plot(t0,np.array([10.,50.]),T,h, rhs)
plot(t0,np.array([15.,30.]),T,h, rhs)
plot(t0,np.array([7.,50.]),T,h, rhs)
plot(t0,np.array([12.,40.]),T,h, rhs)
plot(t0,np.array([70.,60.]),T,h, rhs)
plot(t0,np.array([50.,80.]),T,h, rhs)
plot(t0,np.array([50.,65.]),T,h, rhs)
plot(t0,np.array([30.,90.]),T,h, rhs)

#plt.savefig('B2_A1.png')
plt.show()
