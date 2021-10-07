import numpy as np
import matplotlib.pyplot as plt

# here we define the original RK4 integrator according to the most general RK4 formula
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

# define the function to integrate
def f(x,y):
    return y/(1+x**2)

n = 200; # define number of steps
x = np.linspace(-20,20,n + 1)
y = np.zeros(n + 1)
y[0] = 1 # define initial conditions
h = x[1] - x[0] # define the interval

# implement the original RK4 method.
for i in range(n):
    y[i + 1] = y[i] + rk4_step(f,x[i],y[i],h)

n2 = round(200/11*4)
x2 = np.linspace(-20, 20, n2 + 1)
yy = np.zeros(n2 + 1)
yy[0] = 1
h2 = x2[1] - x2[0]
for i in range(n2):
    yy[i + 1] = rk4_stepd(f,x2[i],yy[i],h2)


n3 = round(200/11*8)
x3 = np.linspace(-20, 20, n3 + 1)
y3 = np.zeros(n3 + 1)
y3[0] = 1
h3 = x3[1] - x3[0]
for i in range(n3):
    y3[i + 1] = rk4_stepd(f,x3[i],y3[i],h3)


y_pred = 4.57605801*np.exp(np.arctan(x))
y_pred2 = 4.57605801*np.exp(np.arctan(x2))
y_pred3 = 4.57605801*np.exp(np.arctan(x3))

plt.clf()
plt.plot(x,y,linewidth = 3,label = 'Original RK4')
plt.plot(x2,yy,'--', label = 'New RK4')
plt.plot(x,y_pred,'-.', label = 'Analytic Solution')
plt.legend()
plt.savefig("RK4.png")
plt.show()

err_h = np.std(y - y_pred)
err_h2 = np.std(yy - y_pred2)
err_h3 = np.std(y3 - y_pred3)
print(f'The error using the normal RK4 method is: {err_h:.3e}')
print(f'The error using the new RK4 method is: {err_h2:.3e}')
print(f'The error using the new RK4 method with 145 steps is: {err_h3:.3e}')


