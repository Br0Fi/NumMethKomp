# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 13:00:36 2014
Edited on Thur May 07 2020

@author: gurevics
edited by Tobias
"""
# Blatt2, SoSe 2020
#--------------------------------------------------------------------
#Solve 1d advection equation with Fourier method
# u_t+c*u_x=0 on x\in[-L/2, L/2] with periodic boundary conditions
# Analytical solution: u(x,t)=u0(x-ct), where
# u0(x,0) is an initial condition.
#--------------------------------------------------------------------
# Movie: avconv -r 25 -i pic%03d.png -b 16M -y filename.mp4

import numpy as np
import matplotlib.pyplot as plt

L = 2.0 * np.pi     #Interval length
N = 128           #Number of points on [-L/2, L/2], N = 2^p
Tend=1200           # t\in[0,Tend]
h=0.01              # Time step
c=0.6               # Propagation velocity

"""x = (np.arange(N)-N/2) * L/N #x-Array		# Interval [-pi,pi]
k = np.zeros(N) #Definition of the wavenumbers		# Alternative definition of k-grid
k[0:N/2+1] = np.arange(N/2+1)*2*np.pi/L
k[N/2+1:] = np.arange(-N/2+1,0)*2*np.pi/L"""
x = (np.arange(N)) * L/N #x-Array 
dx=L/N
k = np.fft.fftfreq(N,dx/(2.*np.pi))

u0= np.exp(- 2*np.pi * (x - L / 4)**2 ) # define an initial condition
uk = np.fft.fft(u0) # uk=fft(u0)

def uexact(L,c,x,h,t): # exact solution of the problem
		return np.exp(-2*np.pi*( np.fmod(np.fmod(x-c*t*h, L)+L,L) -np.pi/2)**2)	# first map on (-2pi,2pi) (because of negative values), then on (0,4pi) and finally [0,2pi)


def rhs(c,k,u): # RHS of the ODE system in the Fourier space
		return -1j*k*c*u

ii=0

plt.ion()  # Force interactive; otherwise each timestep the figure has to be closed manually

for n in range(Tend+1): # RK4 time loop in Fourier space

    k1 =rhs(c,k,uk)
    k2 =rhs(c,k,uk+0.5*h*k1)
    k3 =rhs(c,k,uk+0.5*h*k2)
    k4 =rhs(c,k,uk+h*k3)
    uk = uk + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    if n % 20 == 0: # plot the solution in x-space
        u = np.fft.ifft(uk)
        plt.cla() #clear axis
        plt.plot(x,uexact(L,c,x,h,n),'b') #compare to exact solution
        plt.plot(x,u.real,'r.') # plot the numerical solution
        plt.ylabel("u(x,t)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(n) + ')')
        plt.show() 
        #plt.draw()     #without plt.ion(), use this command instead of plt.show()
        ii=ii+1
        filename = "pic{:03d}.png".format(ii)
        plt.savefig(filename)
        plt.pause(0.1)# pause between frames


plt.ioff()
plt.show() #with this, figure does not close automatically    
