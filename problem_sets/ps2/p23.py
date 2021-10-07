import numpy as np

def mylog2(x):
    # separate out the mantissa and exponent for our x and e
    ma_x, exp_x = np.frexp(x)
    ma_e, exp_e = np.frexp(np.e)

    # evaluate the ma of e and x using our previous log2 function
    ma = np.polynomial.chebyshev.chebval(ma_xj,n)
    me = np.polynomial.chebyshev.chebval(ma_e,n)

    #using a change of basis
    lnx = (ma + exp_x)/(me + exp_e)
    return lnx

# define the interval of evaluate
inv = 100
x = np.linspace(0.5,1,inv)
y = np.log2(x)

x_rescale = np.linspace(-1,1,inv) # rescaling the x to from -1 to 1 to fit the chebyshev polynomial
n = np.polynomial.chebyshev.chebfit(x_rescale,y,7) # define the order of polynomial which gives the uncertainty below 10e-6
che = np.polynomial.chebyshev.chebval(x_rescale,n) #evaluate the log2 at values of rescale

print(np.max(y-che))

x_m= np.linspace(0.1, 10, 1000) # define the interval of evaluation
ln_act = np.log(x_m) # have the actual function
err = abs(mylog2(x_m) - ln_act)
print(err.mean())
