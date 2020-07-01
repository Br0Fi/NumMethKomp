
# command to create movie: avconv -r 25 -i SH-2D%04d.png -b 16M -y SH.mp4

import numpy as np
import matplotlib.pyplot as plt
import os


Lx = 32 #Physical length of the area in x direction
Ly = 32 #Physical length of the area in y direction
Nx = 128 # number of grid points in x direction
Ny = 128 # number of grid points in y direction

x, y = np.meshgrid(np.arange(Nx) * Lx/Nx,np.arange(Ny) * Ly/Ny) #x-Array
kx, ky = np.meshgrid(np.fft.fftfreq(Nx,Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny,Ly/(Ny*2.0*np.pi)))
ksq = kx*kx + ky*ky

# Initial conditions:
csp = 0.0+(np.random.random((Nx,Ny))-0.5)*0.1

eps=0.1
b = -0.1 #0.0
c = 1.0/16

t = 0.0
dt = 0.01
T_End = 500
N_t = int(T_End / dt)

plotEveryNth = 200
#
ck = np.fft.fft2(csp) 

# Linear terms
def rhs_lin(ksq):
    result=eps-(1.0-ksq)**2
    return result

#function for the semiimplicit step
Eu=1.0/(1.0-dt*rhs_lin(ksq))

folders=["5a","5b"]
folder = folders[1]     #select part here!
if not os.path.exists("./"+folder):
    os.makedirs("./"+folder,exist_ok=True)

for i in range(N_t):
    ck=Eu*(ck+dt*(-c*np.fft.fft2(np.fft.ifft2(ck)*np.fft.ifft2(ksq**2*np.fft.fft2(np.fft.ifft2(ck)**2)))-b*np.fft.fft2(np.fft.ifft2(ck)**3)))
    
    if(i % plotEveryNth == 0):
        csp=np.fft.ifft2(ck)
        print("Step {:04d}".format(round(i/plotEveryNth)) )
        plt.suptitle('L='+str(int(Lx))+', N='+str(Nx)+', $\\varepsilon=$'+str(eps)+', b='+str(b))
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
