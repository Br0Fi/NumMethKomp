# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

""" Array init """

x_vec = np.linspace(0,2*np.pi,100)
k_vec = np.fft.fftshift(np.fft.fftfreq(np.size(x_vec), 2*np.pi/np.size(x_vec))) # Frequenzen, nicht omega


""" Funktionen def """

def f(x):
    
    return np.sin(x)

def g(x):
    
    return np.exp(-(x-np.pi)**2)

def FTf(x):
    
    return 0 #außer Delta-Distributionen

def FTg(x):
    
    return np.exp(-x**2/4+1.0j*np.pi*x)/np.sqrt(2)/np.sqrt(2*np.pi)/2*np.pi*np.size(x_vec) #Normierung anpassen

def g_prime(x):
    
    return -2*(x-np.pi)*np.exp(-(x-np.pi)**2)


""" FFT """

FT_f = np.fft.fft(f(x_vec))
FT_g = np.fft.fft(g(x_vec))


""" IFFT """

IFT_f = np.fft.ifft(FT_f)
IFT_g = np.fft.ifft(FT_g)


""" Frequenzen zentrieren """

FT_f = np.fft.fftshift(FT_f)
FT_g = np.fft.fftshift(FT_g)


""" Ableitung durch FFT """

FFTg_prime = np.fft.ifft(1.0j*2*np.pi*k_vec*FT_g)


plt.ion()
plt.rc('text', usetex=True)

""" Vergleich Ursprungsfunktion mit Funktion durch IFT der FT """
plt.figure(1,figsize=(24, 12))
plt.plot(x_vec,f(x_vec),'b',linewidth=2.0)
plt.plot(x_vec,IFT_f,'r',linestyle='--',linewidth=2.0)
plt.plot(x_vec,g(x_vec),'g',linewidth=2.0)
plt.plot(x_vec,IFT_g,'orange',linestyle='--',linewidth=2.0)
plt.ylabel("f(x) bzw. g(x)", fontsize = 24)
plt.xlabel("x", fontsize = 24)
plt.title("Ursprungsfunktions vs. IFT", fontsize = 32)
plt.axis([0, 2*np.pi, np.min(IFT_f), np.max(IFT_f)])
plt.tick_params(labelsize = 24)
plt.legend(["f(x)", "f(x) durch IFFT","g(x)", "g(x) durch IFFT"], fontsize = 24)
plt.savefig("1a.pdf")
plt.show()

""" Vergleich der FFT mit der analytischen Lösung"""

plt.cla()
#plt.plot(k_vec,FTf(k_vec),'b',linewidth=2.0)
plt.vlines([-1,1], 0, 100, 'b')
plt.plot(2*np.pi*k_vec,FT_f,'r',linestyle='--',linewidth=2.0)
plt.plot(k_vec,FTg(k_vec),'g',linewidth=2.0)
plt.plot(2*np.pi*k_vec,FT_g,'orange',linestyle='--',linewidth=2.0)
plt.ylabel("F(x) bzw. G(x)", fontsize = 24)
plt.xlabel("k", fontsize = 24)
plt.title("FFT vs. analytische Lösung", fontsize = 32)
plt.axis([-5, 5, np.min(FTg(x_vec)), np.max(FTg(x_vec))])
plt.tick_params(labelsize = 24)
plt.legend(["F(x) durch FFT", "G(x)","G(x) durch FFT", "F(x)"], fontsize = 12)
plt.savefig("1b.pdf")
plt.show()

""" Vergleich der FFT-Ableitung mit der analytischen Lösung"""

plt.cla()
plt.plot(x_vec,g_prime(x_vec),'b',linewidth=2.0)
plt.plot(x_vec,FFTg_prime,'r',linestyle='--',linewidth=2.0)
plt.ylabel("g\'(x)", fontsize = 24)
plt.xlabel("x", fontsize = 24)
plt.title("Ableitung durch FFT vs. analytische Lösung", fontsize = 32)
plt.axis([0, 2*np.pi, np.min(FFTg_prime), np.max(FFTg_prime)])
plt.tick_params(labelsize = 24)
plt.legend(["g\'(x)","g\'(x) durch FFT"], fontsize = 12)
plt.savefig("1c.pdf")
plt.show()
