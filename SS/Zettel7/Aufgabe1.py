
# command to create movie: avconv -r 25 -i SH-2D%04d.png -b 16M -y SH.mp4

import numpy as np
import matplotlib.pyplot as plt
import os


Lx = 128 #Physical length of the area in x direction
Ly = 128 #Physical length of the area in y direction
Nx = 128 # number of grid points in x direction
Ny = 128 # number of grid points in y direction

x, y = np.meshgrid(np.arange(Nx) * Lx/Nx,np.arange(Ny) * Ly/Ny) #x-Array
kx, ky = np.meshgrid(np.fft.fftfreq(Nx,Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny,Ly/(Ny*2.0*np.pi)))
ksq = kx*kx + ky*ky

# Initial conditions:
csp = 0.0+(np.random.random((Nx,Ny))-0.5)*0.1

kappa = 0.3

t = 0.0
dt = 0.01
T_End = 1000
N_t = int(T_End / dt)

plotEveryNth = 100
#
ck = np.fft.fft2(csp) 

# Linear terms
def rhs(ksq,ck):
    result=-ksq*((kappa*ksq-1)*ck+np.fft.fft2(np.fft.ifft2(ck)**3))
    return result


folders=["kappa1","kappa2","kappa.5","kappa.3","kappa.1"]
folder = folders[3]     #select part here!
if not os.path.exists("./"+folder):
    os.makedirs("./"+folder,exist_ok=True)

for i in range(N_t):
    k1 = rhs(ksq, ck)
    k2 = rhs(ksq, ck + 0.5*dt*k1)    
    k3 = rhs(ksq, ck + 0.5*dt*k2)    
    k4 = rhs(ksq, ck + dt*k3)
    ck = ck + (dt/6.)*(k1 + 2.*k2 + 2.*k3 + k4) 
    
    if(i % plotEveryNth == 0):
        csp=np.fft.ifft2(ck)
        print("Step {:04d}".format(round(i/plotEveryNth)) )
        plt.suptitle('L='+str(int(Lx))+', N='+str(Nx)+', $\\kappa=$'+str(kappa))
        plt.subplot(1,2,1)
        plt.cla()
        plt.imshow(csp.real)
        plt.ylabel("y",fontsize=18)
        plt.xlabel("x",fontsize=18)
        plt.title('$\\psi$(x, y, '+ str(i) + ')')
        plt.axis('off')
        #
        plt.subplot(1,2,2)
        plt.cla()
        plt.imshow(np.abs(np.fft.fftshift(ck)))
        plt.title('Fourier space, $|\\widehat{\\psi}_k|$')
        plt.axis('off')
        #
        
        filename = folder+"/SH-2D{:04d}.png".format(round(i/plotEveryNth))
        plt.savefig(filename)
        plt.close()
