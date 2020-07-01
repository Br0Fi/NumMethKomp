import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import os

Lx = 64.
Ly = 64.
Nx = 256
Ny = 256
T_End = 400
t = 0.0
dt = 0.01
N_t = int(T_End/dt)

x, y = np.meshgrid(np.linspace(0., Lx*(Nx-1.)/Nx, Nx), np.linspace(0., Ly*(Ny-1.)/Ny, Ny))
kx, ky = np.meshgrid(np.fft.fftfreq(Nx,Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny,Ly/(Ny*2.0*np.pi)))
ksq = kx*kx+ky*ky

folder = "A3b"
plotEveryNth = 2000

epsilon = 0.3
delta = 0.0
k = 1.178
rd = 2*np.pi*random.rand() #willkürliches \varphi zwischen 0.0 und 2\pi
c = np.sqrt(epsilon/3)*(np.exp(1j*(rd+k*x))+np.exp(-1j*(rd+k*x)))+(np.random.random((Nx,Ny))-0.5)*0.1
ck = np.fft.fft2(c)

def rhs_lin(ksq):
    return epsilon-(1.0-ksq)**2

if not os.path.exists("./"+folder):
    os.makedirs("./"+folder,exist_ok=True)

for i in range(N_t):
    ck = (ck+dt*(delta*np.fft.fft2(np.fft.ifft2(ck)**2)-np.fft.fft2(np.fft.ifft2(ck)**3)))/(1.-dt*rhs_lin(ksq))
    
    if (i % plotEveryNth == 0):
        c = np.fft.ifft2(ck)
        plt.suptitle("$L=$"+str(int(Lx))+", $N=$"+str(Nx)+", $\\varepsilon=$"+str(epsilon)+", $\\delta=$"+str(delta)+", $k=$"+str(k)+", $\\varphi=$"+str(round(rd, 3)))
        plt.subplot(1,2,1)
        plt.cla()
        plt.imshow(c.real)
        plt.title("$\\psi$(x, y, "+ str(int(i*dt)) + ") in real space:")
        plt.axis("off")
        
        plt.subplot(1,2,2)
        plt.cla()
        plt.imshow(np.abs(np.fft.fftshift(ck)))
        plt.title("$|\\widehat{\\psi}_k|$ in Fourier space:")
        plt.axis("off")
                
        filename = "/A3b-2D{:04d}.png".format(round(i/plotEveryNth))
        plt.savefig(folder+filename)
        plt.close()
