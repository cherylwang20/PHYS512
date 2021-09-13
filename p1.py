import numpy as np

# define the point which we want to evalute
x = 42
# define a function to use


def f1(n):
    return np.exp(n)


def f2(n):
    return np.exp(0.01*n)


# define the interval appropriate for evaluating the integrals, see explanation in doc
eps = 6*10**(-16)
d1 = (27*eps*f1(x)/2/f1(x))**(1/5)
d2 = (27*eps*f2(x)/2/0.01**5/f2(x))**(1/5)


def deriv(n):
    if n == 1:
        x1 = x + 2*d1
        x2 = x + d1
        x3 = x - d1
        x4 = x - 2*d1
        print(f1)
        return (f1(x4) - 8*f1(x3) + 8*f1(x2) - f1(x1))/(12*d1)
    if n == 2:
        x1 = x + 2 * d2
        x2 = x + d2
        x3 = x - d2
        x4 = x - 2 * d2
        return (f2(x4) - 8*f2(x3) + 8*f2(x2) - f2(x1))/(12*d2)


num = np.exp(x)
num2 = 0.01*np.exp(0.01*x)

print(num)

print('The derivative of exp(x) at ', x, 'is' , deriv(1), 'the fractional error is', abs(deriv(1)/num - 1))
print('The derivative of exp(0.01x) at ', x, 'is', deriv(2), 'the error is', abs(deriv(2)/num2 - 1))




