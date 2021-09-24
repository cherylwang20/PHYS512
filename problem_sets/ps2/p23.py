import numpy as np

def mylog2(x):
    ma_x, exp_x = np.frexp(x)
    ma_e, exp_e = np.frexp(np.e)

    ma = np.polynomial.chebyshev.chebval(ma_x,n)
    me = np.polynomial.chebyshev.chebval(ma_e,n)

    lnx = (ma + exp_x)/(me + exp_e)
    return lnx

inv = 100
x = np.linspace(0.5,1,inv)
y = np.log2(x)

# rescaling the x to from -1 to 1 to fit the chebyshev polynomial
x_rescale = np.linspace(-1,1,inv)
n = np.polynomial.chebyshev.chebfit(x_rescale,y,7)
che = np.polynomial.chebyshev.chebval(x_rescale,n)

print(np.max(y-che))


x_m= np.linspace(0.1, 10, 1000)
ln_act = np.log(x_m)
err = abs(mylog2(x_m) - ln_act)
print(err.mean())
