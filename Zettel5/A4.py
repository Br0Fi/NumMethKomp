import numpy as np
import matplotlib.pyplot as plt


# Globale Funktion:
def f(r, x_n):
    return r * (x_n - x_n**2)


# Simulation der Fixpunktbestimmung:
def Spinnennetz(r, x0, n, Dateiname):
    x = np.linspace(0, 1, 1000)
    N = int(n)
    xn = np.zeros(N+1)
    xn[0] = x0
    
    plt.figure(figsize=(10, 7.5))
    plt.plot(x, f(r, x), "-b")
    plt.plot(x, x, "-k")
    plt.plot([x0, x0], [0., x0], "--r") #Verbindungslinie zwischen x-Achse und Gerade bei x0
        
    for i in range(0, N):
        xn[i+1] = f(r, xn[i])
        plt.plot([xn[i], xn[i]], [xn[i], xn[i+1]], "--r") #Vertikale Linien
        plt.plot([xn[i], xn[i+1]], [xn[i+1], xn[i+1]], "--r") #Horizontale Linien
    
    plt.xlabel(r"$x$", fontsize=20)
    plt.ylabel(r"$f(x)$", fontsize=20)
    plt.title(r"Simulation des Systems mit $r={}$".format(r), fontsize=20)
    plt.tick_params(labelsize=16)
    plt.grid(True)
    #plt.savefig(Dateiname + ".png")
    plt.show()
    print("Fixpunkt:", xn[N]) #Zeigt den ermittelten Fixpunkt.
    return


# Methodenaufrufe:
Spinnennetz(1.4, 0.2, 100, "A4.1")
Spinnennetz(2.5, 0.125, 100, "A4.2")
Spinnennetz(2.8, 0.75, 100, "A4.3")
Spinnennetz(3.7, 0.4, 100, "A4.4")


#------------------------------------------------------------------------------


# Berechnet einen Fixpunkt, solange r im Intervall (1,3) liegt:
def Spiderweb(r, x0, n):
    N = int(n)
    xn = np.zeros(N+1)
    xn[0] = x0
    for i in range(0, N):
        xn[i+1] = f(r, xn[i])
    
    return xn[N]


# Plotet die berechneten Fixpunkte in Abhaengigkeit von r:
r = np.linspace(0, 4, 10000)
xS = np.zeros(10000)
for i in range(0, 10000):
    xS[i] = Spiderweb(r[i], 0.9, 100)
plt.figure(figsize=(10, 7.5))
plt.plot(r, xS, ".b")
plt.xlabel(r"$r$", fontsize=20)
plt.ylabel(r"$x_n$", fontsize=20)
plt.title(u"Feigenbaum-Bifurkationsdiagramm", fontsize=20)
plt.tick_params(labelsize=16)
plt.grid(True)
#plt.savefig("Feigenbaum.png")
plt.show()