## RK4-METHOD
## numerical integration of x' = rhs(t,x) over time interval [t0,T] with time step h and initial condition x(t0)=x0
## x can have arbitrary dimension, according to the dimension of x0 and the output of rhs

import numpy as np

def rk4(t0, x0, T, h, rhs):
    N = int((T-t0)/h)
    n = len(x0)

    t = np.linspace(t0,T,N+1)
    x = np.empty((N+1,n))
    x[0] = x0

    for i in range(N):
        k1 = rhs(t[i], x[i])
        k2 = rhs(t[i]+h/2., x[i]+h/2.*k1)
        k3 = rhs(t[i]+h/2., x[i]+h/2.*k2)
        k4 = rhs(t[i]+h, x[i]+h*k3)
        x[i+1] = x[i] + h/6.*(k1 + 2.*(k2 + k3) + k4)
    return t, x
