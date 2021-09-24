import numpy as np
import matplotlib.pyplot as plt

inv = 100

x = np.linspace(0.5,1,inv)
y = np.log2(x)

# rescaling the x to from -1 to 1 to fit the chebyshev polynomial
x_rescale = np.linspace(-1,1,inv)
n = np.polynomial.chebyshev.chebfit(x_rescale,y,8)
che = np.polynomial.chebyshev.chebval(x_rescale,n)

print(np.max(y-che))

def mylog2(x):
    ma_x, exp_x = np.frexp(x)
    ma_e, exp_e = np.frexp(np.e)

    ma = np.polynomial.chebyshev.chebval(ma_x,n)
    me = np.polynomial.chebyshev.chebval(ma_e,n)

    lnx = (ma + exp_x)/(me + exp_e)
    return lnx

x_m= np.linspace(0.1, 10, 20)
ln_act = np.log(10)
err = abs(mylog2(10) - ln_act)
print(err)