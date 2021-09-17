import numpy as np
from scipy import interpolate


# cubic spline interpolation
def lakeshore_c(V, data):

    # extract the initial V and T data from the imported matrix
    V_in = data[:, 1]
    T_in = data[:, 0]

    # we need to inverse the order of both V and T since the cubic spline
    # method only takes in array with increasing values
    spln = interpolate.splrep(V_in[::-1], T_in[::-1])

    # take in the inputted desire V to interpolate the T value(s) we need
    TT_c = interpolate.splev(V, spln)
    return TT_c

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

# bootstrap resampling to evaluate the accuracy of data we interpolated
# we take a portion of the sample as the data for resampling for a single point
# and we find the mean and standard deviation for each data point as a means of
# uncertainty

def error_fun(V,data):
    # we generate a new array to append the new T datasets
    gen_pts = []
    # intake the data set
    V_int = data[:,1]
    T_int = data[:,0]

    # repeat the bootstrap method for a number of times
    for i in range(N_resample):
        # intake the indices of the data sample
        indices = list(range(V_int.size))
        # randomly choose the indices we use to interpolated with a desired N_sample
        to_interp = rng.choice(indices, size = N_sample, replace=False)
        # reverse the order so the indices are in ascending orders
        to_interp = sorted(to_interp, reverse=True)
        # create a new data set points from the chosen new indices
        data_choice = np.vstack((V_int[to_interp],T_int[to_interp])).T
        # input the new T and V data set into our previous interpolation function
        # also input the V data that we want to evaluate the accuracy with
        new_T = lakeshore_c(V, data_choice)
        #append to the array we created
        gen_pts.append(new_T)

    gen_pts = np.array(gen_pts)
    stds = np.std(gen_pts,axis=0) # Calculate std dev at each T accross resamplings
    error_mean = np.mean(stds) # Take mean of std for overall error
    error_std = np.std(stds) # Get std dev of that error.
    return error_std,error_mean
    # print(f"{error_mean = :.3e} +/- {error_std: .3e}")

if __name__ == '__main__':
    # import the data set
    dat = np.loadtxt('lakeshore.txt')
    T = dat[:,0]
    V = dat[:,1]
    # generate a random seed for bootstrap
    rng = np.random.default_rng(seed = 12345)
    # choose the size of resample and sample
    N_resample = len(V)//5
    N_sample = len(V) - N_resample

    # create a list of Vs we want to evaluate through interpolation
    V_i = np.array([0.6,0.8]) # enter value between 0.1 and 1.64
    T_new = lakeshore_c(V_i, dat)
    err = error_fun(V_i, dat) # evaluate the error function

    print('The Temperature at', V_i, 'are',T_new ,'with error', f"{err[1] :.3e} +/- {err[0]: .3e}")

