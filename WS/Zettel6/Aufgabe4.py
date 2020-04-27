import numpy as np
import matplotlib.pyplot as plt

"""Differenzverfahren für DGL der Form x: = p(t)*x. + q(t)*x + r(t)"""

def DV(p, q, r, x_0, x_n, t_0, t_n, N_steps):
    
    delta_t = (t_n - t_0)/N_steps
    
    #Matrix aufstellen
    A = np.zeros((N_steps-1,N_steps-1))
    
    for i in range(N_steps-1):
        
        A[i,i] = -1*(2 + delta_t**2 * q(t_0 + (i+1)*delta_t))
        
        if i < N_steps - 2:
            
            A[i,i+1] = 1 - p(t_0 + (i+1)*delta_t)
            
        if i > 0:
            
            A[i,i-1] = 1 + p(t_0 + (i+1)*delta_t)
            
    
    # Inhomogenität (Vektor) aufstellen
    b = np.zeros((N_steps-1))
    
    for i in range(N_steps-1):
        
        b[i] = delta_t**2 * r(t_0 + (i+1)*delta_t)
        
        if i == N_steps - 2:
            
            b[i] -= x_n * (1 - delta_t/2 * p(t_0 + (i+1)*delta_t))
            
        if i == 0:
            
            b[i] -= x_0 * (1 + delta_t/2 * p(t_0 + (i+1)*delta_t))
    
    return np.linalg.solve(A,b)


"""Eingabeerte und -funktionen"""
    
N_steps = 100
    
t_0 = -1
t_n = 1
x_0 = 0
x_n = 0
    
def p(t):
        
    return 0
    
def q(t):
        
    return -(1+t**2)
    
def r(t):
        
    return -1
    
x = np.concatenate(([x_0], DV(p, q, r, x_0, x_n, t_0, t_n, N_steps), [x_n]))
    
t = np.linspace(t_0, t_n, N_steps + 1)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t,x)

plt.title(r'RWP über Differenzverfahren', fontsize=14, pad=15)
plt.xlabel('$t$', fontsize=12)
plt.ylabel('$x$', fontsize=12)

plt.savefig('Aufgabe4.pdf', bbox_inches='tight')
plt.show()