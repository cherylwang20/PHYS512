import numpy as np

def ndiff(fun,x,full = False):
    dx = np.zeros(len(x))
    deriv = np.zeros(len(x))
    error = np.zeros((len(x)))
    error2 = np.zeros(len(x))
    for i in range(len(x)):
        eps = 7 * 10 ** (-17)  # what should be the optimized eps???
        dx[i] = (eps) ** (1 / 3) * x[i]
        deriv[i] = (fun(x[i] + dx[i]) - fun(x[i] - dx[i]))/(2*dx[i])
        error[i] = eps*fun(x[i])/dx[i] + dx[i]**2/6*fun(x[i]) #how to express the truncation error term?
        error2[i]= eps*fun(x[i])/dx[i]

    if full == False:
        return deriv
    if full == True:
        return deriv,dx,abs(error),abs(error2)

fun = np.exp
x = [1]


print('The numerical derivatives of',fun,' at', x, 'is/are', ndiff(fun,x,False))

# what is the exact error term??