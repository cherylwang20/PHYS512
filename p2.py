import numpy as np


eps = 7 * 10 ** (-17)
def ndiff(fun,x,full = False):
    dx = np.zeros(len(x))
    deriv = np.zeros(len(x))
    error = np.zeros((len(x)))
    error2 = np.zeros(len(x))
    for i in range(len(x)):
        dx[i] = d(fun,x[i])[0]
        deriv[i] = (fun(x[i] + dx[i]) - fun(x[i] - dx[i]))/(2*dx[i])
        error[i] = abs(d(fun,x[i])[1])

    if full == False:
        return deriv
    if full == True:
        return deriv,dx,abs(error)

def d(f, x):
    # define the fifth derivative function
    h = eps ** (1 / 4)
    func = (f(x+2*h) - 3*f(x+h) + 3*f(x) - f(x- h))/ h ** 3
    delt = (3 * eps * abs(f(x)) / abs(func)) ** (1 / 3)
    err = delt**2*abs(func)/6 + eps*abs(f(x))/delt
    return delt , err

#define the function and x values that we want to evaluate
fun = np.exp
x = [0.3,1.4,3.2]

def f2(n):
    return 1/(n**2+1)

#er = abs(ndiff(fun,x,True)[0] /np.exp(x) - 1)
#print(er, ndiff(fun,x,True)[2])
print('The numerical derivatives of',fun,' at', x, 'is/are', ndiff(fun,x,True)[0],' with intervals', ndiff(fun,x,True)[1],'and uncertainty', ndiff(fun,x,True)[2] )
print('The numerical derivatives of',f2,' at', x, 'is/are', ndiff(f2,x,True)[0])