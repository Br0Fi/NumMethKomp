
# command to create movie: avconv -r 25 -i GL-2D%04d.png -b 16M -y GL2D.mp4

import numpy as np
import matplotlib.pyplot as plt


Lx = 32. #Physical length of the area in x direction
Ly = 32. #Physical length of the area in x direction
Nx = 256  # number of grid points in x direction
Ny = 256  # number of grid points in y direction
T_End = 200
tt=50*2
dt=0.05
N_t = int(T_End / dt)
print(N_t)
#
delta=1.0
epsilon=0.1
#
x, y = np.meshgrid(np.arange(Nx) * Lx/Nx,np.arange(Ny) * Ly/Ny) #x-array
kx, ky = np.meshgrid(np.fft.fftfreq(Nx,Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny,Ly/(Ny*2.0*np.pi)))
ksq = kx*kx + ky*ky

A0= 1.0*0 + (np.random.random((Nx,Ny))-0.5)*0.01#0.1 *2
A=A0
Ahat=np.fft.fft2(A)
#
#
q=epsilon-(-ksq+1)**2
coef1=((1+q*dt)*np.exp(q*dt)-1-2*dt*q)/(dt*q**2)
coef2=(-np.exp(q*dt)+1+dt*q)/(dt*q**2)
#
for n in range(N_t):
    Nn=np.fft.fft2(delta*np.fft.ifft2(Ahat)**2-np.fft.ifft2(Ahat)**3)
    if n==0:
        Nn1=Nn
    Ahat=Ahat*np.exp(q*dt)+Nn*coef1+Nn1*coef2
    Nn1=Nn
    A=np.fft.ifft2(Ahat)
    if (n % tt == 0):
        print("Step {:04d}".format(round(n/tt)))
        fig1=plt.figure(1, figsize=(16, 12))
        fig1.suptitle('Frosen States in 2D CGLE:$\\,\\delta$='+str(delta)+','+'$\\,\\epsilon=$'+str(epsilon), fontsize=20)
        plt.subplot(1,2,1)
        plt.cla()
#        im1 = plt.imshow(A.real, vmin=-1., vmax=1.)
        im1 = plt.imshow(A.real)
        plt.colorbar(im1)
        plt.ylabel("y",fontsize=18)
        plt.xlabel("x",fontsize=18)
        plt.title('Re(A(x, y, '+ str(n) + '))')
        #plt.axis('off')
        #
        plt.subplot(1,2,2)
        plt.cla()
        im2 = plt.imshow(np.abs(A))
        plt.colorbar(im2)
        plt.title('|A|')
        plt.ylabel("y",fontsize=18)
        plt.xlabel("x",fontsize=18)
#        plt.axis('off')
        #
        filename = "GL-2D{:04d}.png".format(round(n/tt))
        plt.savefig(filename)
        plt.close()
        
        