import numpy as np
from scipy import interpolate


# polynomial interpolation
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

# cubic spline interpolation
def lakeshore_c(V, data):
    V_in = data[:, 1]
    T_in = data[:, 0] # for later comparsion

    spln = interpolate.splrep(V_in[::-1], T_in[::-1])
    TT_c = interpolate.splev(V, spln)
    return TT_c


#bootstrap resampling

def error_fun(V,data):
    gen_pts = []
    V_int = data[:,1]
    T_int = data[:,0]
    for i in range(100):
        indices = list(range(V_int.size))
        to_interp = rng.choice(indices, size = N_sample, replace=False)
        to_interp = sorted(to_interp, reverse=True)
        data_choice = np.vstack((V_int[to_interp],T_int[to_interp])).T
        #to_check = [i for i in indices if not (i in to_interp)]
        new_T = lakeshore_c(V, data_choice)
        gen_pts.append(new_T)

    gen_pts = np.array(gen_pts)
    stds = np.std(gen_pts,axis=0)
    error_mean = np.mean(stds)
    error_std = np.std(stds)
    return error_std,error_mean
    # print(f"{error_mean = :.3e} +/- {error_std: .3e}")

if __name__ == '__main__':
    dat = np.loadtxt('lakeshore.txt')
    T = dat[:,0]
    V = dat[:,1]
    rng = np.random.default_rng(seed = 12345)
    N_resample = len(V)//5
    N_sample = len(V) - N_resample
    V_i = np.array([0.6,0.8]) #enter value between 0.1 and 1.64
    print(lakeshore_c(V_i, dat))
    err = error_fun(V_i,dat)
    print('The Temperature at', V_i, 'are',lakeshore(V_i, dat) ,'with error', f"{err[1] :.3e} +/- {err[0]: .3e}")

