import numpy as np
import matplotlib.pyplot as plt

alpha = 1.5

bins=np.linspace(0,10,501)
cents=0.5*(bins[1:]+bins[:-1])

mylor=1/(1+cents**2)*3
myexp = np.exp(-alpha*cents)
mygau = np.exp(-alpha*cents**2)*3
mypow = cents**(-alpha)


plt.plot(cents, myexp, label = 'Exponential')
plt.plot(cents, mylor, label = 'Lorenztian')
plt.plot(cents, mygau, label = 'Gaussian')
plt.plot(cents - 1, mypow, label = 'Power')
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

eff_loz = len(t_use)/len(t)
print(f'The efficiency of the Lorenztian bounding distribution is: {eff_loz}.')

aa,bb=np.histogram(t_use,bins)
b = 0.5*(bb[1:] + bb[:-1])
aa = aa/aa.sum()
exp = myexp/myexp.sum()
plt.bar(b, aa, 0.05, label = 'Exponential Deviates')
plt.plot(b, exp, 'r', label = 'Predicted')
plt.title('Rejection Method via Lorenz')
plt.legend()
plt.savefig('Reject_Lor.png',dpi = 300, bbox_inches = 'tight')
plt.show()

# rejection method via power laws

q=np.random.rand(n)
t_pow=(q)**(1/(1-alpha))
y_pow = t_pow**(-alpha)*np.random.rand(n)

accept_pow=y_pow<np.exp(-alpha*t_pow)
t_use_p=t_pow[accept_pow]
print(t_use_p)

aa_pow,bb=np.histogram(t_use_p,bins + 1)
aa_pow = aa_pow/aa_pow.sum()
plt.bar(cents , aa_pow, 0.05, label = 'Exponential Deviates')
plt.plot(cents , exp, 'r', label = 'Predicted')
plt.title('Rejection Method via Power Law')
plt.legend()
plt.savefig('Reject_Pow.png',dpi = 300, bbox_inches = 'tight')
plt.show()

eff_pow = len(t_use_p)/len(t_pow)
print(f'The efficiency of the Power Law bounding distribution is: {eff_pow}.')