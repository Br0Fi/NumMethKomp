## RKF(4,5)-METHOD
## numerical integration of x' = rhs(t,x) over time interval [t0,T] with time step h and initial condition x(t0)=x0
## x can have arbitrary dimension, according to the dimension of x0 and the output of rhs

import numpy as np

# alpha is just the left column of beta
beta = [[0.], [1/4, 1/4], [3/8, 3/32, 9/32], [12/13, 1932/2197, -7200/2197, 7296/2197], [1., 439/216, -8, 3680/513, -845/4104],
        [1/2, -8/27, 2, -3544/2565, 1859/4104, -11/40]]
cestar = [0., 16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
ce = [0., 25/216, 0, 1408/2565, 2197/4104, -1/5, 0] # for simplicities sake c_0 = c*_0 = 0 is added

# calculate one step of RK4 or RK5
def calc_x(t, x, h, rhs, stage): # stage = 4 oder 5
    k = np.empty([len(beta), len(x)]) # e.g. k_1 -> k_0
    xp = np.empty(len(x))
    for i1 in range(len(k)):
        #set t value:
        tmod = t+beta[i1][0]*h
        #set x values for each k using the last k:
        xmod = np.zeros(len(x))
        for i2 in range(1, i1):
            xmod += beta[i1][i2] * k[i1-1]
        # xn = xn + h* (sum_i=1^s-1 (beta[i1][i]*k_i1-1))
        xmod = x + xmod * h # nope, think for 1 and higher i1
        k[i1] = rhs(tmod, xmod)

    xp = np.zeros(len(x))
    for i in range(len(k[0])):
        if stage==4:
            xp += ce[i] * k[i]
        elif stage==5:
            xp += cestar[i] * k[i]
        else:
            print("calc_x error: stage must be 4 or 5, but is:")
            print(stage)
    xp = x + xp*h
    return xp

#test = np.array([[1,2,3],[4,5,6],[7,8,9]])
#test[0][1] = 42
#print(test[0])

# RUNGE KUTTA FEHLBERG
# N is passed instead of T so that we can work with a fixed length array
# it could be implemented to have eps_0 be an array with a value for each variable
# you could use one epsilon_0 for all variables and as long as one error is larger than eps_0
# step size is increased, but this makes the epsilon_0/epsilon for h_new be weird.
# so here just the "main" variable x[0] is checked
def rkf(t0, x0, N, h_start, rhs, eps_0, beta_ctrl):
    n = len(x0)
    t = np.empty(N)
    x = np.empty((N,n))
    x[0] = x0
    h = h_start
    t[0] = t0
    for i in range(N-1):
        # calculate the next step
        x1 = calc_x(t[i], x[i], h, rhs, 4) # for rk4
        x2 = calc_x(t[i], x[i], h, rhs, 5) # for rk5
        t[i] = t0 + h
        eps = np.abs(x1[0]-x2[0]) # see above the def
        if eps<eps_0: # check if any eps is larger than eps_0
            if(eps!=0.):
                h = beta_ctrl * h * (eps_0/eps)**(1/6) # p+1 = 6, weil p=5 I guess?
            else:
                h = h_start
            h = min(h, h_start) # otherwise step-size gets infinite, but this doesn't seem great either.
            t[i+1] = t[i] + h
            x[i+1] = x2
        else:
            h = beta_ctrl * h * (eps_0/eps)**(1/5) # see above
            i -= 1 # redo the calculation in this case.
    return t, x

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
