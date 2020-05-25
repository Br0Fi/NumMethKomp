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
Tend= 20000          # t\in[0,Tend]
h=0.000001              # Time step
c1=5
c2=4

"""x = (np.arange(N)-N/2) * L/N #x-Array		# Interval [-pi,pi]
k = np.zeros(N) #Definition of the wavenumbers		# Alternative definition of k-grid
k[0:N/2+1] = np.arange(N/2+1)*2*np.pi/L
k[N/2+1:] = np.arange(-N/2+1,0)*2*np.pi/L"""
x = np.arange(N) * L/N - L/2 #x-Array from -pi to pi
dx=L/N
k = np.fft.fftfreq(N,dx/(2.*np.pi))
u0= c1**2/2/np.cosh(c1*(x+2)/2)**2 + c2**2/2/np.cosh(c2*(x+1)/2)**2
# define an initial condition

uk = np.fft.fft(u0) # uk=fft(u0)
wk = uk

def dealiase(u): #how is that spelled?
	return u

def rhs(k,wk,tt): # RHS in the Fourier space
	expk3 = np.exp(-1j*k**3*tt)
	return -3*1j*k* expk3* np.fft.fft(np.fft.ifft(expk3*wk)**2)

def rk4f(k,uk,tt):
    k1 =rhs(k,uk,tt)
    k2 =rhs(k,uk+0.5*h*k1,tt)
    k3 =rhs(k,uk+0.5*h*k2,tt)
    k4 =rhs(k,uk+h*k3,tt)
    return uk + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

ii=0

plt.ion()  # Force interactive; otherwise each timestep the figure has to be closed manually

for n in range(Tend+1): # RK4 time loop in Fourier space

    tt = ii*h
    wk = rk4f(k,wk,tt)

    if n % 50 == 0: # plot the solution in x-space

        u = np.fft.ifft(np.exp(1j * k**3 * tt)*wk)
        plt.cla() #clear axis
        #plt.ylim(-1,1)
        plt.plot(x,u.real,'r.') # plot the numerical solution
        plt.ylabel("u(x,t)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(tt) + ')')
        plt.show()
        ii=ii+1
        filename = "pic{:03d}.png".format(ii)
        #plt.savefig("pics/"+filename)
        plt.pause(0.01)# pause between frames


plt.ioff()
#plt.show() #with this, figure does not close automatically
