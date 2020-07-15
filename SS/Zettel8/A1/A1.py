import numpy as np
import matplotlib.pyplot as plt

L = 120
N = 128
du = 5.
dv = 12.
#a = 2.
#a = 2.5
a = 3.
#b = 9.
#b = 12.
#b = 15.
b = 18.
dt = 0.01
T_End = 150
N_t = int(T_End / dt)
plotEveryNth = 250

x, y = np.meshgrid(np.linspace(0., L*(N-1.)/N, N), np.linspace(0., L*(N-1.)/N, N))
kx, ky = np.meshgrid(np.fft.fftfreq(N,L/(N*2.0*np.pi)), np.fft.fftfreq(N,L/(N*2.0*np.pi)))
ksq = kx*kx+ky*ky

u = a + (np.random.random((N,N))-0.5)*0.1
v = b/a + (np.random.random((N,N))-0.5)*0.1
uk = np.fft.fft2(u)
vk = np.fft.fft2(v)

for i in range(N_t):
    uk = (uk+dt*np.fft.fft2(np.fft.ifft2(uk)**2*np.fft.ifft2(vk)+a))/(1.+dt*(du*ksq+b+1.))
    vk = (vk+dt*b*uk-dt*np.fft.fft2(np.fft.ifft2(uk)**2*np.fft.ifft2(vk)))/(1.+dt*dv*ksq)
    u = np.fft.ifft2(uk)
    v = np.fft.ifft2(vk)
    
    if(i % plotEveryNth == 0):
        print("Step %04d" %(i/plotEveryNth))
        plt.suptitle("L="+str(int(L))+", N="+str(N)+", a="+str(a)+", b="+str(b))
        plt.subplot(1,2,1)
        plt.cla()
        plt.imshow(u.real)
        plt.ylabel("y",fontsize=18)
        plt.xlabel("x",fontsize=18)
        plt.title("u(x, y, "+ str(dt*i) + ")")
        plt.axis("off")
        
        plt.subplot(1,2,2)
        plt.cla()
        plt.imshow(v.real)
        plt.title("v(x, y, "+ str(dt*i) + ")")
        plt.axis('off')
        
        filename = "RD-2D%04d.png" % (i/plotEveryNth)
        plt.savefig(filename)
        plt.close()