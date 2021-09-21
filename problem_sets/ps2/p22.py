import numpy as np

def lorentz(x):
    return 1/(1+x**2)

def integrate_adaptive(fun,a,b,tol,extra = None):
    x = np.linspace(a,b,5)
    y = fun(x)
    dx = (a - b)/(len(x)-1)
    return
