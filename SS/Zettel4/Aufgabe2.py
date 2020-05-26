# Movie: avconv -r 25 -i pic%03d.png -b 16M -y filename.mp4

import numpy as np
import matplotlib.pyplot as plt

Lx = 2.0 * np.pi     #Interval length x
Ly = 2.0 * np.pi     #Interval length y
Nx = 128           #Number of points on [-L/2, L/2], N = 2^p x
Ny = 128           #Number of points on [-L/2, L/2], N = 2^p y
Tend=1200           # t\in[0,Tend]
h=0.01              # Time step
cx=0.6               # Propagation velocity x
cy=0.6               # Propagation velocity y

"""x = (np.arange(N)-N/2) * L/N #x-Array		# Interval [-pi,pi]
k = np.zeros(N) #Definition of the wavenumbers		# Alternative definition of k-grid
k[0:N/2+1] = np.arange(N/2+1)*2*np.pi/L
k[N/2+1:] = np.arange(-N/2+1,0)*2*np.pi/L"""
x , y = np.meshgrid(np.arange(Nx) * Lx/Nx, np.arange (Ny) * Ly/Ny)
dx=Lx/Nx
dy=Ly/Ny
kx, ky = np.meshgrid(np.fft.fftfreq(Nx, Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny, Ly/(Ny*2.0*np.pi)))

u0= np.exp(-2*np.pi*( ( x -Lx / 4) **2+(y - Ly / 4) **2) )# define an initial condition
uk = np.fft.fft(u0) # uk=fft(u0)


def rhs(cx,kx,cy,ky,u): # RHS of the ODE system in the Fourier space
		return -1j *( kx *cx+ky*cy ) *u

ii=0

plt.ion()  # Force interactive; otherwise each timestep the figure has to be closed manually

for n in range(Tend+1): # RK4 time loop in Fourier space

    k1 =rhs(cx,kx,cy,ky,uk)
    k2 =rhs(cx,kx,cy,ky,uk+0.5*h*k1)
    k3 =rhs(cx,kx,cy,ky,uk+0.5*h*k2)
    k4 =rhs(cx,kx,cy,ky,uk+h*k3)
    uk = uk + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    if n % 20 == 0: # plot the solution in x-space
        u = np.fft.ifft2(uk)
        plt.cla() #clear axis
        plt.imshow(np.abs(u), cmap='coolwarm')
        plt.colorbar()
        plt.ylabel("y")
        plt.xlabel("x")
        plt.title('u(x,y,'+ str(n) + ')')
        #plt.draw()     #without plt.ion(), use this command instead of plt.show()
        ii=ii+1
        filename = "pic{:03d}.png".format(ii)
        plt.savefig(filename)
        plt.show() 
        plt.pause(0.1)# pause between frames


plt.ioff()
plt.show() #with this, figure does not close automatically    