import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft, ifft

N = 1024
k = 10.4
x = np.arange(N)
y = np.sin(2*np.pi*k*x/N)
y_ft = np.abs(np.fft.fft(y))

fn = np.arange(0, 1, 1/N)


ke  = x
J = np.complex(0, 1)
f1 = [0]*N
f2 = [0]*N
for i in range(len(ke)):
    f1[i] = (1 - np.exp(-2*np.pi*J*(ke[i] + k)))/(1-np.exp(-2*np.pi*J*(ke[i] + k)/(N)))
    f2[i] = (1 - np.exp(-2 * np.pi * J * (ke[i] - k))) / (1 - np.exp(-2 * np.pi * J * (ke[i] - k) /(N)))
f = np.abs((np.array(f1) - np.array(f2))/(2*J))
#f = f/np.sum(f)

plt.plot(fn, y_ft,label = 'True')
plt.plot(fn, f,label = 'Analytic')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.legend()
plt.savefig(f'Analyic_Sine_k = {k}_N = {N}.png',dpi = 300, bbox_inches = 'tight')
plt.show()

res = np.std(np.abs(f - y_ft))

print(f'The residue between the functions are: {res}, with N = {N}.' )
#the closer to the end, the larger the residue.
plt.plot(fn, f - y_ft)
plt.title('Residual between the two functions')
plt.savefig(f'Residual.png',dpi = 300, bbox_inches = 'tight')
plt.show()

win = 0.5 - 0.5*np.cos(2*np.pi*x/N)
y_win = y*win
y_win_ft = np.abs(np.fft.fft(y_win))


plt.plot(fn, y_win_ft, label = 'Windowed Function')
plt.legend()
plt.savefig(f'Sine_Window.png',dpi = 300, bbox_inches = 'tight')
plt.show()
