import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


# here we use the legender method as our own integrator

def get_legendre_weight(n):
    x = np.linspace(-1,1,n+1)
    P = np.polynomial.legendre.legvander(x,n)
    Pinv = np.linalg.inv(P)
    coeffs = Pinv[0,:]
    return coeffs*n

coeffs = get_legendre_weight(50)

def legender_int(a,x):
    my_int = np.empty(len(a))
    for i in range(len(a)):
        # define the function to integrate
        y = (a[i] - R*x)/(R**2 + a[i]**2 - 2*R*a[i]*x)**(3/2)
        dx = x[1] - x[0]
        my_int[i]  = np.sum(coeffs*y)*dx
    return my_int

x = np.linspace(-1,1,len(coeffs))
# define the constants used in the integration
R = 1.4
a = np.linspace(0.01*R,5*R,100)
sig = 100 #*1.6*10**(-19) #eV
eps = 8.8 #*10**(-12)

int_leg = legender_int(a,x)*sig*R**2/2/eps

# second integration method using scipy.quad
intt = np.empty(len(a))
for i in range(len(a)):
    y2 = lambda x: (a[i] - R*x)/(R**2 + a[i]**2 - 2*R*a[i]*x)**(3/2)
    intt[i] = integrate.quad(y2,-1,1)[0]*sig*R**2/2/eps

print(sig*R**2/2/eps)
plt.figure(1)
plt.ylim(-20,60)
plt.xlim(0,max(a))
plt.text(0.8*R,40,'R',color = 'red')
plt.text(1.8*R,1.2*sig*R**2/2/eps,'Sigma*R^2/2/eps_0')
plt.plot(a,int_leg,label = "Legendre Integral")
plt.plot(a, intt,'-.', label = "Scipy.quad Integral")
plt.axvline(x = R,color = 'red',linestyle = 'dashed')
plt.axhline(y = sig*R**2/2/eps, color = 'green',linestyle = '-.')
plt.legend()
plt.show()


