import numpy as np

def lorentz(x):
    return 1/(1+x**2)


n = 0

def integrate_adaptive(fun, a, b, tol, extra = None):
    x = np.linspace(a,b,5)
    y = fun(x)
    dx = (b - a) / (len(x) - 1)
    global n
    if extra is None:
        area1 = 2 * dx * (y[0] + 4 * y[2] + y[4]) / 3
        area2 = dx * (y[0] + 4 * y[1] + 2 * y[2] + 4 * y[3] + y[4]) / 3
        err = np.abs(area1 - area2)
    else:
        area1 = 2*dx*(extra[0] + 4*fun(x[2]) + extra[2])/3
        area2 = dx * (extra[0] + 4*fun(x[1]) + 2 * fun(x[2]) + 4 * fun(x[3]) + extra[2])/3
        err = np.abs(area1 - area2)
        n = n + 4
    if err < tol:
        return area2
    else:
        left = integrate_adaptive(fun,a,x[2],tol/2, y)
        right = integrate_adaptive(fun,x[2],b,tol/2, y)
    return left+right

a = -100
b = 100
ans = integrate_adaptive(lorentz,a,b,1e-7)
act = np.arctan(b) - np.arctan(a)
print(n)
print(ans)
print(abs(ans-act))

