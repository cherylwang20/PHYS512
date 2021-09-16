import numpy as np

# define functions to use
def f1(n):
    return np.exp(n)

def f2(n):
    return np.exp(0.01 * n)

def deriv(d, f, x):
    x1 = x + 2 * d(f, x)
    x2 = x + d(f, x)
    x3 = x - d(f, x)
    x4 = x - 2 * d(f, x)
    return (f(x4) - 8 * f(x3) + 8 * f(x2) - f(x1)) / (12 * d(f, x))

x = 2
# define the interval appropriate for evaluating the integrals
eps = 6 * 10 ** (-16)
def d(f, x):
    # define the fifth derivative function
    h = eps ** (1 / 6)
    func = (f(x+4*h) - 5*f(x+3*h) + 10*f(x+2*h) - 10*f(x+h) + 5*f(x) - f(x- h))/ h ** 5
    return (27 * eps * f(x) / (2 * abs(func))) ** (1 / 5)


# this function estimates the derivative, taking in the desire interval d and function f, at point x
est = deriv(d, f1, x)
num = np.exp(x)
num2 = 0.01 * np.exp(0.01 * x)

print(est, num)

print('The derivative of exp(x) at ', x, 'is', deriv(d, f1, x), ',the fractional error is',
      abs(deriv(d, f1, x) / num - 1))
print('The derivative of exp(0.01x) at ', x, 'is', deriv(d, f2, x), ',the fractional error is',
      abs(deriv(d, f2, x) / num2 - 1))
