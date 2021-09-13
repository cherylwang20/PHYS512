import numpy as np
from matplotlib import pyplot as plt

dat = np.loadtxt('lakeshore.txt')


def lakeshore(V, data):
    T_in = data[:, 0]
    V_in = data[:, 1]

    ind = np.argmin(V_in > V) # find the argument of the smallest number greater than V
    print(ind)
    V_use = V_in[ind - 1:ind + 3]
    T_use = T_in[ind - 1:ind + 3]
    p = np.polyfit(V_use, T_use, 3)
    T_out = np.polyval(p, V)
    return T_out


print(lakeshore(0.135480, dat))
# make it into an array &error term
# use