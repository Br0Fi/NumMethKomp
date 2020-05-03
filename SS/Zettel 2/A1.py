import numpy as np
import matplotlib.pyplot as plt

def u_initial(x):
    return np.exp(-2.*np.pi*(x-np.pi/2.)**2)

def f(k, c, uk):
    return -1.*1j*k*c*uk

N = 200
time = 1020
h = 0.01
c = 0.6
x = np.linspace(0, 2*np.pi, N)
k = np.fft.fftfreq(N, 1./N)
u0 = u_initial(x)
uk = np.fft.fft(u0)

for i in range(0, time+1):
    k1 = f(k, c, uk)
    k2 = f(k, c, uk+h*k1/2.)
    k3 = f(k, c, uk+h*k2/2.)
    k4 = f(k, c, uk+h*k3)
    uk = uk+(h/6.)*(k1+2.*(k2+k3)+k4)
    
    if i%20==0:
        u = np.fft.ifft(uk)
        plt.figure(figsize=(10, 7.5))
        plt.plot(x, u.real, ".r", zorder=3, label="Numerische Loesung")
        plt.plot(x, u_initial(np.fmod(np.fmod(x-h*c*i, 2.*np.pi)+2*np.pi, 2*np.pi)), zorder=2, label="Analytische Loesung")
        plt.legend(loc="best", fontsize=18)
        plt.xlabel(r"$x$", fontsize=20)
        plt.xlim(0., 2.*np.pi)
        plt.ylabel(r"$u(x,t=$"+str(round(h*i,1))+"$)$", fontsize=20)
        plt.ylim(-0.1, 1.3)
        plt.tick_params(labelsize=16)
        plt.grid(True)
        title = r"Loesung mit $N=$%3d, $h=$%4.2f und $c=$%3.1f" % (N, h, c)
        plt.title(title, fontsize=20)
        filename = "plot%4d.png" % i
        plt.savefig(filename)
        plt.show()