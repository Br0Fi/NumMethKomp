import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

fig_size = (7, 7)

def mask(row, rule):
    # stack vertically shifted versions of array to get
    # the pattern for each cell
    y = np.vstack((np.roll(row, 1), row,
                   np.roll(row, -1))).astype(np.int8)
    # give a representation from 0 to 7 to each pattern:
    #  0,   1,   2,  ...,  6,   7
    # 111, 110, 101, ..., 001, 000
    u = np.array([[4], [2], [1]])
    z = np.sum(y * u, axis=0).astype(np.int8)
    # return array for next step given by rule array
    return rule[7 - z]

def calc_field(a, rule):
    for i in range(0, len(a)-1):
        a[i+1] = mask(a[i], rule)
    return a

#parameters and initial condition
T = 101
size = 101
rule30 = np.array([0, 0, 0, 1, 1, 1, 1, 0])
rule90 = np.array([0, 1, 0, 1, 1, 0, 1, 0])
    
a = np.zeros((T, size)) # space on one axis, time on the other
a[0][int(size/2)] = 1

calc_field(a, rule30)
plt.figure(figsize=fig_size)
plt.imshow(a,interpolation='nearest',cmap=plt.cm.gray_r)
plt.show()

a = np.zeros((T, size)) # space on one axis, time on the other
a[0][int(size/2)] = 1

calc_field(a, rule90)
plt.figure(figsize=fig_size)
plt.imshow(a,interpolation='nearest',cmap=plt.cm.gray_r)
plt.show()
