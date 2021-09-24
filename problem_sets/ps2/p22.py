import numpy as np

# define the functio of our choice
def lorentz(x):
    return 1/(1+x**2)


n = 0 # count the # of function calls we save

def integrate_adaptive(fun, a, b, tol, extra = None):
    x = np.linspace(a,b,5)
    y = fun(x)
    dx = (b - a) / (len(x) - 1)
    global n
    if extra is None: # first time calling
        area1 = 2 * dx * (y[0] + 4 * y[2] + y[4]) / 3
        area2 = dx * (y[0] + 4 * y[1] + 2 * y[2] + 4 * y[3] + y[4]) / 3
        err = np.abs(area1 - area2) # compare the two methods of integration
    else:
        area1 = 2*dx*(extra[0] + 4*fun(x[2]) + extra[2])/3 # we reuse two of the points from previous function
        area2 = dx * (extra[0] + 4*fun(x[1]) + 2 * fun(x[2]) + 4 * fun(x[3]) + extra[2])/3 # reuse two points fro previous function
        err = np.abs(area1 - area2)
        n = n + 4 # increase the count of function calls we save
    if err < tol: # if the error is smaller than tolerance, return the area
        return area2
    else:
        left = integrate_adaptive(fun,a,x[2],tol/2, y) #we store the function we already evaluted into the next round
        right = integrate_adaptive(fun,x[2],b,tol/2, y)
    return left+right

# define the interval we want to evalute
a = -100
b = 100
ans = integrate_adaptive(lorentz,a,b,1e-7) # define the tolerance level
act = np.arctan(b) - np.arctan(a) # compare it to the actual function
print(n)
print(ans)
print(abs(ans-act))

