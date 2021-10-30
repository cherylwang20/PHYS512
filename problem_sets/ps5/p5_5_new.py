import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft, ifft

k = 40.555
fs = 100
N = 256
n = np.arange(0, N)
x = n/fs

y = np.sin(2*np.pi*k*x/N)

NFFT = 1024
nVAL = np.arange(0,NFFT)
yft = np.abs(fft(y,NFFT))
yft = yft/np.sum(yft)
fre = x*fs/N

ke = fre*N
J = np.complex(0, 1)
f1 = [0]*N
f2 = [0]*N
for i in range(len(ke)):
    f1[i] = (1 - np.exp(2*np.pi*J*(k - ke[i])))/(1-np.exp(2*np.pi*J*(k - ke[i])/(N)))
    f2[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] + k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] + k) /(N)))
f_1= np.abs((np.array(f1) - np.array(f2))/(2*J))
f = f_1/np.sum(f_1)



plt.plot(nVAL, np.fft.fftshift(yft),label = 'true')
#plt.plot(fre,f,label = 'analytic')
plt.legend()
plt.title(f'Non-Integer Sine Function with k = {k}')
#plt.savefig(f'Analyic_Sine_{k}.png',dpi = 300, bbox_inches = 'tight')
plt.show()

#
# f1_w = [0]*N
# f2_w = [0]*N
# for i in range(len(ke)):
#     f1_w[i] = (1 - np.exp(4*np.pi*J*(k - ke[i])))/(1-np.exp(4*np.pi*J*(k - ke[i])/(N)))
#     f2_w[i] = (1 - np.exp(-4 * np.pi * J * (ke[i] + k))) / (1 - np.exp(-4 * np.pi * J * (ke[i] + k) /(N)))
# f_w = np.abs((np.array(f1) - np.array(f2))/(2*J))/4
#
#
# f_final = f_1/2 + f_w/4
# f_final = f_final/np.sum(f_final)
# plt.plot(fre,yft,label = 'true')
# plt.plot(fre,np.fft.fftshift(f_final),label = 'analytic')
# plt.legend()
# plt.title(f'Non-Integer Sine Function with N = {N} data points')
# plt.savefig(f'Analyic_Sine_{N}_win.png',dpi = 300, bbox_inches = 'tight')
# plt.show()
#
# w = 0.5 - 0.5*np.cos(2*np.pi*x/N)
# w_ft = np.fft.fft(w)
# print(w_ft[0:2])
# plt.plot(fre,w_ft)
# plt.show()