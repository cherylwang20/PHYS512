import numpy as np
from matplotlib import pyplot as plt

k = 1/3
x = np.arange(0,9,0.03)
N = len(x)
y = np.sin(2*np.pi*x*k/N)
yft = np.abs(np.fft.fftshift(np.fft.fft(y)))
yft = yft/np.sum(yft)

ke  = np.arange(0,N)
J = np.complex(0, 1)
f1 = [0]*N
f2 = [0]*N
for i in range(len(ke)):
    f1[i] = (1 - np.exp(-2*np.pi*J*(ke[i] + k)))/(1-np.exp(-2*np.pi*J*(ke[i] + k)/(N+1)))
    f2[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] - k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] - k) /(N+1)))
f = np.abs((np.array(f1) - np.array(f2))/(2*J))
f = f/np.sum(f)

# plt.plot(x,yft,label = 'true')
# plt.plot(x,np.fft.fftshift(f),label = 'analytic')
# plt.legend()
# plt.xlim(4,5)
# plt.title(f'Non-Integer Sine Function with N = {N} data points')
# plt.savefig(f'Analyic_Sine_{N}.png',dpi = 300, bbox_inches = 'tight')
# plt.show()


win=0.5-0.5*np.cos(2*x*np.pi/N)
f1_w = [0]*N
f2_w = [0]*N
for i in range(len(ke)):
    f1_w[i] = (1 - np.exp(-2*np.pi*J*(ke[i] + k)))/(1-np.exp(-2*np.pi*J*(ke[i] + k)/(N+1)))
    f2_w[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] - k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] - k) /(N+1)))
f = np.abs((np.array(f1) - np.array(f2))/(2*J))
f = f/np.sum(f)

