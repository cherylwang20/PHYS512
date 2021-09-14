import numpy as np

dat = np.loadtxt('lakeshore.txt')


def lakeshore(V, data):
    T_in = data[:, 0]
    V_in = data[:, 1]


    if isinstance(V,float) == True:
        ind = np.argmin(V_in > V)  # find the argument of the smallest number greater than V
        V_use = V_in[ind - 1:ind + 3]
        T_use = T_in[ind - 1:ind + 3]
        p = np.polyfit(V_use, T_use, 3)
        T_out = np.polyval(p, V)
    else:
        T_out = np.empty(len(V))
        for i in range(len(V)):
            ind = np.argmin(V_in > V[i])  # find the argument of the smallest number greater than V
            V_use = V_in[ind - 1:ind + 3]
            T_use = T_in[ind - 1:ind + 3]
            p = np.polyfit(V_use, T_use, 3)
            T_out[i] = np.polyval(p, V[i])
    return T_out

V_i = [0.58,0.130,1.45,1.3,1.25666]

print('The Temperature at', V_i, 'are',lakeshore(V_i, dat),'correspondingly.')
# What is the error term?
