import numpy as np
import matplotlib.pyplot as plt

def rk4_step(fun,x,y,h):
    k0 = fun(x,y)*h
    k1 = fun(x + h/2, y + k0/2)*h
    k2 = fun(x + h/2, y + k1/2)*h
    k3 = fun(x + h, y + k2)*h
    return (k0 + 2*k1 + 2*k2 + k3)/6

def rk4_stepd(fun,x,y,h):
        st1_h = y + rk4_step(fun,x,y,h)
        st2_h = y + rk4_step(fun,x,y,h/2)
        st3_h = st2_h + rk4_step(fun,x+h/2, st2_h ,h/2)
        err = (st1_h - st3_h)/15
        return st3_h + err

def f(x,y):
    return y/(1+x**2)

n = 200
x = np.linspace(-20,20,n + 1)
y = np.zeros(n + 1)
y[0] = 1
h = x[1] - x[0]

for i in range(n):
    y[i + 1] = y[i] + rk4_step(f,x[i],y[i],h)

n2 = round(200/11*4)
x2 = np.linspace(-20, 20, n2 + 1)
yy = np.zeros(n2 + 1)
yy[0] = 1
h2 = x2[1] - x2[0]
for i in range(n2):
    yy[i + 1] = rk4_stepd(f,x2[i],yy[i],h2)


y_pred = 4.57605801*np.exp(np.arctan(x))
y_pred2 = 4.57605801*np.exp(np.arctan(x2))

plt.clf()
plt.plot(x,y)
plt.plot(x2,yy,'--')
plt.plot(x,y_pred,'-.')
plt.show()

err_h = np.std(y - y_pred)
err_h2 = np.std(yy - y_pred2)
print(err_h,err_h2)


