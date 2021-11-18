import numpy as np
import matplotlib.pyplot as plt

n=100000
u=np.random.rand(n)
v=(np.random.rand(n)*2-1)*0.8

r=v/u

alpha = 1.5

accept=u<np.sqrt(np.exp(-alpha*r))
t=r[accept]

bins=np.linspace(0,10,501)
cents=0.5*(bins[1:]+bins[:-1])

myexp = np.exp(-alpha*cents)

aa,bb=np.histogram(t,bins)
b = 0.5*(bb[1:] + bb[:-1])
aa = aa/aa.sum()
exp = myexp/myexp.sum()
plt.bar(b, aa, 0.05, label = 'Exponential Deviates')
plt.plot(cents, exp, 'r', label = 'Predicted')
plt.title('Rejection Method via Ratio-of-Uniform Generator ')
plt.legend()
plt.savefig('Reject_ru.png',dpi = 300, bbox_inches = 'tight')
plt.show()

eff_ru = len(t)/len(r)
print(f'The efficiency of Ratio-of-Uniform Generator is {eff_ru}')