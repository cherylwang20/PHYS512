import numpy as np
import matplotlib.pyplot as plt

def rk4_step(fun,x,y,h):
    k0 = fun(x,y)*h
    k1 = fun(x + h/2, y + k0/2)*h
    k2 = fun(x + h/2, y + k1/2)*h
    k3 = fun(x + h, y + k2)*h
    return (k0 + 2*k1 + 2*k2 + k3)/6

def f(x,y):
    return y/(1+x**2)

n = 200
x = np.linspace(-20,20,n + 1)
y = np.zeros(n + 1)
y[0] = 1
h = x[1] - x[0]

for i in range(n):
    y[i + 1] = y[i] + rk4_step(f,x[i],y[i],h)

y_pred = 4.57605801*np.exp(np.arctan(x))

plt.clf()
plt.plot(x,y)
plt.plot(x,y_pred,'-.')
plt.show()

err_h = np.std(y - y_pred)
print(err_h)

def rk4_stepd(fun,x,y,h):

    return