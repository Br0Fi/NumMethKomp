# Blatt3, SoSe 2020
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
N = 256           #Number of points on [-L/2, L/2], N = 2^p
Tend=1500           # t\in[0,Tend]
h=0.01              # Time step
nus=[0.01,0.007,0.004,0.001]
nu=nus[1]

"""x = (np.arange(N)-N/2) * L/N #x-Array		# Interval [-pi,pi]
k = np.zeros(N) #Definition of the wavenumbers		# Alternative definition of k-grid
k[0:N/2+1] = np.arange(N/2+1)*2*np.pi/L
k[N/2+1:] = np.arange(-N/2+1,0)*2*np.pi/L"""
x = np.arange(N) * 2*L/N - L #x-Array from -2pi to 2pi
dx=L/N
k = np.fft.fftfreq(N,dx/(2.*np.pi))
u0= np.sin(x) # define an initial condition, mind mistake in exercise sheet

uk = np.fft.fft(u0) # uk=fft(u0)

def dealiase(u): #how is that spelled?
	return u

def rhs(nu,k,uk): # RHS in the Fourier space
	uk = dealiase(uk)
	dduk = 1j*k*uk # calculate derivative
	NN = -1*np.fft.ifft(uk) * np.fft.ifft(dduk) # calculate the nonlinearity
	NNk = dealiase(np.fft.fft(NN)) # transform and dealiase
	Lk = -1 * nu * k*k * uk # I think there was an error in the lecture at this step (nu^2)
	return Lk - NNk # signs irritating

def rk4f(nu,k,uk):
    k1 =rhs(nu,k,uk)
    k2 =rhs(nu,k,uk+0.5*h*k1)
    k3 =rhs(nu,k,uk+0.5*h*k2)
    k4 =rhs(nu,k,uk+h*k3)
    return uk + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

ii=0

plt.ion()  # Force interactive; otherwise each timestep the figure has to be closed manually

for n in range(Tend+1): # RK4 time loop in Fourier space

    uk = rk4f(nu,k,uk)

    if n % 20 == 0: # plot the solution in x-space

        u = np.fft.ifft(uk)
        plt.cla() #clear axis
        plt.ylim(-1,1)
        plt.plot(x,u.real,'r.') # plot the numerical solution
        plt.ylabel("u(x,t)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(n) + '),  ' + "nu =" + str(nu))
        plt.show()
        #plt.draw()     #without plt.ion(), use this command instead of plt.show()
        ii=ii+1
        filename = "pic{:03d}.png".format(ii)
        #plt.savefig("pics/"+filename)
        plt.pause(0.1)# pause between frames


plt.ioff()
plt.show() #with this, figure does not close automatically
