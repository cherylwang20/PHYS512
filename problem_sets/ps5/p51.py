import numpy as np
from matplotlib import pyplot as plt

def con_shift(array,x):
    return
x = np.arange(-10,10,.05)
y = np.exp(-0.5*x**2/(1.5**2))
shift = len(x)/2
kvec = np.arange(len(x))
J = np.complex(0,1)
print(len(kvec))
yft_shift = np.fft.fft(y)*np.exp(-2*np.pi*J*kvec*shift/len(x))
y_shift = np.real(np.fft.ifft(yft_shift))

plt.plot(x,y)
plt.plot(x,y_shift)
plt.savefig("gauss_shift.png",dpi = 300, bbox_inches = 'tight')
plt.show()
