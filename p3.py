import numpy as np
from scipy import interpolate

dat = np.loadtxt('lakeshore.txt')


def lakeshore(V, data):
    T_in = data[:, 0] # for later comparsion
    V_in = data[:, 1]

    if isinstance(V,float) == True:
        ind = np.argmin(V_in > V)  # find the argument of the smallest number greater than V
        if ind < 1: # set boundary conditions
            V_use = V_in[ind:ind + 4]
            T_use = T_in[ind:ind + 4]
            p = np.polyfit(V_use, T_use, 3)
            T_out = np.polyval(p, V)
        else:
            V_use = V_in[ind - 1:ind + 3]
            T_use = T_in[ind - 1:ind + 3]
            p = np.polyfit(V_use, T_use, 3)
            T_out = np.polyval(p, V)
    else:
        T_out = np.empty(len(V))
        for i in range(len(V)):
            ind = np.argmin(V_in > V[i])  # find the argument of the smallest number greater than V
            if ind < 1:  # set boundary conditions
                V_use = V_in[ind:ind + 4]
                T_use = T_in[ind:ind + 4]
                p = np.polyfit(V_use, T_use, 3)
                T_out = np.polyval(p, V)
            else:
                V_use = V_in[ind - 1:ind + 3]
                T_use = T_in[ind - 1:ind + 3]
                p = np.polyfit(V_use, T_use, 3)
                T_out[i] = np.polyval(p, V[i])
    return T_out



def lakeshore_c(V, data):
    V_in = data[:, 1]
    T_in = data[:, 0] # for later comparsion

    spln = interpolate.splrep(V_in[::-1], T_in[::-1])
    TT_c = interpolate.splev(V, spln)    
    return TT_c

def ndiff(fun,x):
    dx = np.zeros(len(x))
    deriv = np.zeros(len(x))
    error = np.zeros((len(x)))
    error2 = np.zeros(len(x))
    for i in range(len(x)):
        eps = 7 * 10 ** (-17)  # what should be the optimized eps???
        dx[i] = (eps) ** (1 / 3) * x[i]
        deriv[i] = (fun(x[i] + dx[i],dat) - fun(x[i] - dx[i],dat))/(2*dx[i])
    return deriv

   
V_i = [0.6,0.8]
print(lakeshore_c(V_i, dat))
ndiff(lakeshore_c,V_i)

print('The Temperature at', V_i, 'are',lakeshore(V_i, dat))
# What is the error term? use np.std over all the points known.
