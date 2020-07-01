import numpy as np
import matplotlib.pyplot as plt
import os

Lx = 32.
Ly = 32.
Nx = 256
Ny = 256
T_End = 300
t = 0.0
dt = 0.01
N_t = int(T_End/dt)

x, y = np.meshgrid(np.linspace(0., Lx*(Nx-1.)/Nx, Nx), np.linspace(0., Ly*(Ny-1.)/Ny, Ny))
kx, ky = np.meshgrid(np.fft.fftfreq(Nx,Lx/(Nx*2.0*np.pi)), np.fft.fftfreq(Ny,Ly/(Ny*2.0*np.pi)))
ksq = kx*kx+ky*ky

folder = "A4a"
plotEveryNth = 1500

epsilon = 0.3
delta = 0.0
c = 0.0+(np.random.random((Nx,Ny))-0.5)*0.1
ck = np.fft.fft2(c)
px = .5*Lx
py = .5*Ly
ex = 15.
#ey = ex
ey = 7.5
epsilon_bar = np.array([[0.] * Nx for i in range(Ny)])
for i in range(Ny):
    for j in range(Nx):
        if (((j*Lx/Nx-px)/ex)**2+((i*Ly/Ny-py)/ey)**2 > 1.):
            epsilon_bar[i][j] = -3.*epsilon

if not os.path.exists("./"+folder):
    os.makedirs("./"+folder,exist_ok=True)

for i in range(N_t):
    ck = (ck+dt*np.fft.fft2(epsilon_bar*np.fft.ifft2(ck)+delta*np.fft.ifft2(ck)**2-np.fft.ifft2(ck)**3))/(1.-dt*(epsilon-(1.0-ksq)**2))
    
    if (i % plotEveryNth == 0):
        c = np.fft.ifft2(ck)
        plt.suptitle("$L=$"+str(int(Lx))+", $N=$"+str(Nx)+", $\\varepsilon=$"+str(epsilon)+", $\\delta=$"+str(delta))
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
                
        #filename = "/A4a-2D-Kreis{:04d}.png".format(round(i/plotEveryNth))
        filename = "/A4a-2D-Ellipse{:04d}.png".format(round(i/plotEveryNth))
        plt.savefig(folder+filename)
        plt.close()