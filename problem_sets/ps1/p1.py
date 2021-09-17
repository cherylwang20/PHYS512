import numpy as np

# define functions to test, we have exp(x) and exp(0.01x)
def f1(n):
    return np.exp(n)

def f2(n):
    return np.exp(0.01 * n)

# this function uses the four point method to acquire the derivative at the
# desire point x using function f and the optimal interval d
# the four point method uses points +/- h and +/- 2h around x
# we take the order to the fifth order error and factor out
# the second, third, and fourth derivative.
def deriv(d, f, x):
    x1 = x + 2 * d(f, x)
    x2 = x + d(f, x)
    x3 = x - d(f, x)
    x4 = x - 2 * d(f, x)
    return (f(x4) - 8 * f(x3) + 8 * f(x2) - f(x1)) / (12 * d(f, x))

# this function returns the optimal interval of h to take to minimize the truncation
# and round-off error. taking the zeros of the derivatives of the sum of truncation and
# round-off error of h, we have the value of h for the minimum error.
def d(f, x):
    # define the fifth derivative function
    h = eps ** (1 / 6)
    func = (f(x+4*h) - 5*f(x+3*h) + 10*f(x+2*h) - 10*f(x+h) + 5*f(x) - f(x- h))/ h ** 5
    return (27 * eps * f(x) / (2 * abs(func))) ** (1 / 5)

if __name__ == '__main__':
    # define the point we want to evaluate
    x = 2
    # define the error interval for x + h
    eps = 6 * 10 ** (-16)

    est1 = deriv(d, f1, x)
    est2 = deriv(d, f2, x)

    # manaully check the error of the derivatives
    num = np.exp(x)
    num2 = 0.01 * np.exp(0.01 * x)

    print('The derivative of exp(x) at ', x, 'is', est1, ',the fractional error is',
          abs(deriv(d, f1, x) / num - 1))
    print('The derivative of exp(0.01x) at ', x, 'is',est2 , ',the fractional error is',
          abs(deriv(d, f2, x) / num2 - 1))
