import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft, ifft

k = 1/3
fs = 1000
nCyl = 5
oversampleRate = 30

x = np.arange(0,nCyl/k - 1/fs, 1/fs )
y = np.sin(2*np.pi*k*x)

yft = np.abs(np.fft.fftshift(np.fft.fft(y)))
yft = yft/np.sum(yft)


ke  = np.arange(0,fs)
J = np.complex(0, 1)
f1 = [0]*fs
f2 = [0]*fs
for i in range(len(ke)):
    f1[i] = (1 - np.exp(-2*np.pi*J*(ke[i] + k)))/(1-np.exp(-2*np.pi*J*(ke[i] + k)/(fs+1)))
    f2[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] - k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] - k) /(fs+1)))
f = np.abs((np.array(f1) - np.array(f2))/(2*J))
f = f/np.sum(f)



fre = np.arange(0, len(x))/len(x)
plt.plot(fre, yft,label = 'true')
#plt.plot(np.fft.fftshift(f),label = 'analytic')
# plt.plot(np.arange(0,N), yftt)
# plt.plot(np.arange(0,N), f)
# plt.legend()
# plt.xlim(-3000,5000)
# plt.title(f'Non-Integer Sine Function with N = {N} data points')
# #plt.savefig(f'Analyic_Sine_{N}.png',dpi = 300, bbox_inches = 'tight')
plt.show()
#
#
# win=0.5-0.5*np.cos(2*x*np.pi/N)
# f1_w = [0]*N
# f2_w = [0]*N
# for i in range(len(ke)):
#     f1_w[i] = (1 - np.exp(-2*np.pi*J*(ke[i] + k)))/(1-np.exp(-2*np.pi*J*(ke[i] + k)/(N+1)))*win
#     f2_w[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] - k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] - k) /(N+1)))*win
# f_win = np.abs((np.array(f1) - np.array(f2))/(2*J))
# f_win = f_win/np.sum(f_win)
#
# plt.plot(x,yft,label = 'true')
# plt.plot(x,np.fft.fftshift(f_win),label = 'analytic')
# plt.legend()
# plt.xlim(4,5)
# plt.title(f'Non-Integer Sine Function with N = {N} data points')
# plt.savefig(f'Analyic_Sine_{N}_win.png',dpi = 300, bbox_inches = 'tight')
# plt.show()