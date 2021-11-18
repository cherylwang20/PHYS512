import numpy as np
import matplotlib.pyplot as plt

alpha = 0.5

bins=np.linspace(0,10,501)
cents=0.5*(bins[1:]+bins[:-1])

mylor=1/(1+cents**2)*3
myexp = np.exp(-alpha*cents)
mygau = np.exp(-alpha*cents**2)
mypow = cents**(-alpha)


plt.plot(cents, myexp, label = 'Exponential')
plt.plot(cents, mylor, label = 'Lorenztian')
plt.plot(cents, mygau, label = 'Gaussian')
plt.plot(cents, mypow, label = 'Power')
plt.ylim(0, 4)
plt.title('Bounding Distribution')
plt.legend()
plt.savefig('Bounding Distribution.png',dpi = 300, bbox_inches = 'tight')
plt.show()

# rejetion method via Lorenztian
def lorentzians(n):
    q=np.pi*(np.random.rand(n)-0.5)
    return np.tan(q)

n=10000000
t=lorentzians(n)
y=1.5/(1+t**2)*np.random.rand(n)*2

accept=y<np.exp(-alpha*t)
t_use=t[accept]

aa,bb=np.histogram(t_use,bins)
aa = aa/aa.sum()
exp = myexp/myexp.sum()
plt.plot(cents, aa, '*', label = 'Exponential Deviates')
plt.plot(cents, exp, 'r', label = 'Predicted')
plt.title('Rejection Method via Lorenz')
plt.legend()
plt.savefig('Reject_Lor.png',dpi = 300, bbox_inches = 'tight')
plt.show()
