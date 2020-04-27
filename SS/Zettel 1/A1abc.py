import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
# A1a:
# Functions:
def f(x):
    return np.sin(x)

def g(x):
    return np.exp(-1. * (x - np.pi)**2)

def A1a(parameter, function, function_name):
    first, last, N = parameter
    
    # Fourier magic:
    x = np.linspace(first, last, int(N))
    signal = function(x)
    k = np.fft.fftfreq(int(N), 1/N)
    fourier = np.fft.fft(signal)
    signal_reconstructed = np.fft.ifft(fourier)
    
    # Plotting:
    plt.figure(figsize=(10, 7.5))
    plt.plot(x, signal)
    plt.xlabel("$x$", fontsize=20)
    plt.xlim(first, last)
    plt.ylabel(function_name+"$(x)$", fontsize=20)
    plt.title("The original signal "+function_name+"$(x)$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)
    plt.savefig("A1a_original_signal_"+function_name[1]+".png")
    plt.show()
    
    plt.figure(figsize=(10, 7.5))
    plt.scatter(k, fourier)
    plt.xlabel("$k$", fontsize=20)
    plt.ylabel(r"$\mathcal{F}[$"+function_name+"$(x)]$", fontsize=20)
    plt.title("The frequency spectrum of "+function_name+"$(x)$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)
    plt.savefig("A1a_frequency_spectrum_"+function_name[1]+".png")
    plt.show()
    
    plt.figure(figsize=(10, 7.5))
    plt.plot(x, signal_reconstructed)
    plt.xlabel("$x$", fontsize=20)
    plt.xlim(first, last)
    plt.ylabel(r"$\mathcal{F}^{-1}[\mathcal{F}[$"+function_name+"$(x)]]$", fontsize=20)
    plt.title("The reconstructed signal "+function_name+"$(x)$", fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)
    plt.savefig("A1a_reconstructed_signal_"+function_name[1]+".png")
    plt.show()



# Using the functions to create the plots:
parameter = [0., 2*np.pi, 128.]
A1a(parameter, f, "$f$")
A1a(parameter, g, "$g$")


#------------------------------------------------------------------------------
# A1c:
# Calculating the derivative of g(x) by using the Fourier transform:
x = np.linspace(0., 2*np.pi, 128)
g_values = g(x)
k = np.fft.fftfreq(128, 1./128.)
fourier_g_values = np.fft.fft(g_values)
tmp = np.zeros(128, dtype=np.complex128)
for i in range(0, 128):
    tmp[i] = 1j * k[i] * fourier_g_values[i]
derivative_g_values = np.fft.ifft(tmp)

# Plotting:
plt.figure(figsize=(10, 7.5))
plt.plot(x, derivative_g_values)
plt.xlabel("$x$", fontsize=20)
plt.xlim(0., 2*np.pi)
plt.ylabel(r"$\frac{d}{dx}g(x)$", fontsize=20)
plt.title("The derivative of $g(x)$", fontsize=20)
plt.tick_params(labelsize=16)
plt.grid(True)
plt.savefig("A1c_derivative_g.png")
plt.show()