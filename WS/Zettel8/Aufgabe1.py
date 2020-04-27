# Blatt 8; Aufgabe 1

import matplotlib.pyplot as plt
import numpy as np
import os


Nt = 1800
dt = 0.05

dx = 0.1
L = 20. 
Nx = np.int(2*L/dx+1)

c = 0.2 # Velocity
alp = dt / dx # Courant number
print("The Courant number alpha = {:.4g}".format(alp))

x = (np.arange(Nx)-(Nx-1)/2) * dx #x-Array [-L L]
t = (np.arange(Nt)) * dt #t-Array [0 T]

#print(x)

u=np.zeros(shape=(Nt+1, Nx)) # solution array

#Anfangsbedingungen für verschiedene Aufgabenteile

def fA(x):  # initial condition
    return 4*np.arctan(np.exp(x/np.sqrt(1-c**2)))

def gA(x):  # initial condition
    return -2*c/np.sqrt(1-c**2)/np.cosh(x/np.sqrt(1-c**2))

def fB(x):  # initial condition
    return 4*np.arctan(np.exp(-x/np.sqrt(1-c**2)))

def gB(x):  # initial condition
    return -2*c/np.sqrt(1-c**2)/np.cosh(x/np.sqrt(1-c**2))

def fC(x):  # initial condition
    return 4*np.arctan(np.exp((x+L/2)/np.sqrt(1-c**2)))+4*np.arctan(np.exp((x-L/2)/np.sqrt(1-c**2)))

def gC(x):  # initial condition
    return -2*c/np.sqrt(1-c**2)/np.cosh((x+L/2)/np.sqrt(1-c**2))+2*c/np.sqrt(1-c**2)/np.cosh((x-L/2)/np.sqrt(1-c**2))

def fD(x):  # initial condition
    return 4*np.arctan(np.exp((x+L/2)/np.sqrt(1-c**2)))+4*np.arctan(np.exp(-(x-L/2)/np.sqrt(1-c**2)))

def gD(x):  # initial condition
    return -2*c/np.sqrt(1-c**2)/np.cosh((x+L/2)/np.sqrt(1-c**2))-2*c/np.sqrt(1-c**2)/np.cosh((x-L/2)/np.sqrt(1-c**2))

def fE(x):  # initial condition
    return 0

def gE(x):  # initial condition
    return 4*np.sqrt(1-c**2)/np.cosh(x*np.sqrt(1-c**2))

#print(u[0,0])

#
A=(1-alp**2)*np.diag(np.ones(Nx))+alp**2/2*np.roll(np.diag(np.ones(Nx)),-1,axis=1)+alp**2/2*np.roll(np.diag(np.ones(Nx)),1,axis=1) 

A[0,1] += alp**2/2
A[Nx-1,Nx-2] += alp**2/2
A[0,Nx-1]=0
A[Nx-1,0]=0


#Aufgabenteil (a)

u[0,:]=fA(x)
gam = gA(x)


if not os.path.exists("./picsA"):
    os.makedirs("./picsA")

plt.ion()
plt.rc('text', usetex=True) # use LaTeX notation

for j in range (Nt):
    
    bet = np.sin(u[j,:])
    
    if (j > 0):
        u[j+1,:]=2*np.dot(A,u[j,:])-u[j-1,:]-dt**2*bet # calculate the next time step
    elif (j == 0):
        u[j+1,:]=np.dot(A,u[j,:])+dt*gam-dt**2/2*bet

    if j % 60 == 0:
        plt.figure(1,figsize=(24, 12))
        plt.cla()
        plt.plot(x,u[j,:],'b',linewidth=2.0)
        plt.ylabel("u(x)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(j) + ')')
        plt.axis([-L, L, np.min(u[j,:]), np.max(u[j,:])])
        filename = "./picsA/pic{:04d}.png".format(j)
        plt.savefig(filename)
        plt.pause(0.02)

plt.ioff()
plt.show()




#Aufgabenteil (b)

u[0,:]=fB(x)
gam = gB(x)


if not os.path.exists("./picsB"):
    os.makedirs("./picsB")

plt.ion()
plt.rc('text', usetex=True) # use LaTeX notation

for j in range (Nt):
    
    bet = np.sin(u[j,:])
    
    if (j > 0):
        u[j+1,:]=2*np.dot(A,u[j,:])-u[j-1,:]-dt**2*bet # calculate the next time step
    elif (j == 0):
        u[j+1,:]=np.dot(A,u[j,:])+dt*gam-dt**2/2*bet

    if j % 60 == 0:
        plt.figure(1,figsize=(24, 12))
        plt.cla()
        plt.plot(x,u[j,:],'b',linewidth=2.0)
        plt.ylabel("u(x)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(j) + ')')
        plt.axis([-L, L, np.min(u[j,:]), np.max(u[j,:])])
        filename = "./picsB/pic{:04d}.png".format(j)
        plt.savefig(filename)
        plt.pause(0.02)

plt.ioff()
plt.show()


#Aufgabenteil (c)

u[0,:]=fC(x)
gam = gC(x)


if not os.path.exists("./picsC"):
    os.makedirs("./picsC")

plt.ion()
plt.rc('text', usetex=True) # use LaTeX notation

for j in range (Nt):
    
    bet = np.sin(u[j,:])
    
    if (j > 0):
        u[j+1,:]=2*np.dot(A,u[j,:])-u[j-1,:]-dt**2*bet # calculate the next time step
    elif (j == 0):
        u[j+1,:]=np.dot(A,u[j,:])+dt*gam-dt**2/2*bet

    if j % 60 == 0:
        plt.figure(1,figsize=(24, 12))
        plt.cla()
        plt.plot(x,u[j,:],'b',linewidth=2.0)
        plt.ylabel("u(x)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(j) + ')')
        plt.axis([-L, L, np.min(u[j,:]), np.max(u[j,:])])
        filename = "./picsC/pic{:04d}.png".format(j)
        plt.savefig(filename)
        plt.pause(0.02)

plt.ioff()
plt.show()




#Aufgabenteil (d)

u[0,:]=fD(x)
gam = gD(x)


if not os.path.exists("./picsD"):
    os.makedirs("./picsD")

plt.ion()
plt.rc('text', usetex=True) # use LaTeX notation

for j in range (Nt):
    
    bet = np.sin(u[j,:])
    
    if (j > 0):
        u[j+1,:]=2*np.dot(A,u[j,:])-u[j-1,:]-dt**2*bet # calculate the next time step
    elif (j == 0):
        u[j+1,:]=np.dot(A,u[j,:])+dt*gam-dt**2/2*bet

    if j % 60 == 0:
        plt.figure(1,figsize=(24, 12))
        plt.cla()
        plt.plot(x,u[j,:],'b',linewidth=2.0)
        plt.ylabel("u(x)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(j) + ')')
        plt.axis([-L, L, np.min(u[j,:]), np.max(u[j,:])])
        filename = "./picsD/pic{:04d}.png".format(j)
        plt.savefig(filename)
        plt.pause(0.02)

plt.ioff()
plt.show()



#Aufgabenteil (e)

u[0,:]=fE(x)
gam = gE(x)


if not os.path.exists("./picsE"):
    os.makedirs("./picsE")

plt.ion()
plt.rc('text', usetex=True) # use LaTeX notation

for j in range (Nt):
    
    bet = np.sin(u[j,:])
    
    if (j > 0):
        u[j+1,:]=2*np.dot(A,u[j,:])-u[j-1,:]-dt**2*bet # calculate the next time step
    elif (j == 0):
        u[j+1,:]=np.dot(A,u[j,:])+dt*gam-dt**2/2*bet

    if j % 60 == 0:
        plt.figure(1,figsize=(24, 12))
        plt.cla()
        plt.plot(x,u[j,:],'b',linewidth=2.0)
        plt.ylabel("u(x)")
        plt.xlabel("x")
        plt.title('u(x,'+ str(j) + ')')
        plt.axis([-L, L, np.min(u[j,:]), np.max(u[j,:])])
        filename = "./picsE/pic{:04d}.png".format(j)
        plt.savefig(filename)
        plt.pause(0.02)

plt.ioff()
plt.show()