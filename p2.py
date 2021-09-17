import numpy as np

# we first write a simple function which uses the two point method
# we create a string of arrays for the derivatives, the optimal intervals
# and the error of each evalutions
def ndiff(fun,x,full = False):
    dx = np.zeros(len(x))
    deriv = np.zeros(len(x))
    error = np.zeros((len(x)))
    error2 = np.zeros(len(x))

    #with the range of x, we evaluate the derivative
    for i in range(len(x)):
        # we call another function to evaluate the optimal dx to use for each x
        dx[i] = d(fun,x[i])[0]
        deriv[i] = (fun(x[i] + dx[i]) - fun(x[i] - dx[i]))/(2*dx[i])
        error[i] = abs(d(fun,x[i])[1])

    if full == False:
        return deriv
    if full == True:
        return deriv,dx,abs(error)

# we write a function to evaluate the optimal dx for each point
# and the estimated error of truncation and round-off
# a detailed derivation could be found in the attached document
def d(f, x):
    # choose the optimal h with trials
    h = eps ** (1 / 4)
    # derive the function for the third derivative use to estimated dx
    func = (f(x+2*h) - 3*f(x+h) + 3*f(x) - f(x- h))/ h ** 3
    # evaluate the optimal delta value for each input x and function f
    delt = (3 * eps * abs(f(x)) / abs(func)) ** (1 / 3)
    # calculate the round-off and truncation error
    err = delt**2*abs(func)/6 + eps*abs(f(x))/delt
    return delt , err

# define a lorentz function for testing
def f2(n):
    return 1/(n**2+1)


if __name__ == '__main__':
    eps = 7 * 10 ** (-17)

    # define the x array and function we want to evaluate
    fun = np.exp
    x = [0.3,1.4,3.2]

    print('The numerical derivatives of',fun,' at', x, 'is/are', ndiff(fun,x,True)[0],' with intervals', ndiff(fun,x,True)[1],'and uncertainty', ndiff(fun,x,True)[2] )
    print('The numerical derivatives of the lorenz function at', x, 'is/are', ndiff(f2,x,True)[0])