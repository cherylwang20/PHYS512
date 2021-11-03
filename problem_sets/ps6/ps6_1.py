import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import json

exec(open("simple_read_ligo.py").read())
f = open(r'C:\Users\wangc\PHYS512\problem_sets\ps6\LOSC_Event_tutorial\BBH_events_v3.json',)

# define the window function that is flat in the middle
def window(n,m):
    x = np.linspace(-np.pi, np.pi, m)
    win = 0.5 + 0.5*np.cos(x)
    mm = m//2
    win_flat = np.ones(n-m)
    win_flat = np.concatenate((win[:mm],win_flat, win[-mm:]))
    return win_flat

data = json.load(f)

print(json.dumps(data, indent = 4, sort_keys=True))

H1, H2, H3, H4 = strain_H


n  = len(strain_H[0])
win = window(len(strain_H[0]), n//5)
th_ft = np.abs(np.fft.rfft(th[0]))
tl_ft = np.abs(np.fft.rfft(tl[0]))
strain_ft = np.fft.rfft(strain_H*win)
sft = np.sum(np.abs(strain_ft),axis = 0)/4

print(sft)

for i in range(4):
    plt.loglog(np.abs(strain_ft[i]))
plt.show()

# Noise Model?
Nft_1 = np.abs(sft)**2
Nft = Nft_1.copy()
# smooth out the noise model:
for i in range(10):
    Nft = (Nft + np.roll(Nft, 1) + np.roll(Nft, -1))/3

plt.loglog(Nft_1, label = 'unsmoothed')
plt.loglog(Nft, label = 'smoothed')
plt.legend()
plt.show()


sft_white = strain_ft/np.sqrt(Nft)
th_white_ft = np.fft.rfft(th*win)/np.sqrt(Nft)
th_white = np.fft.irfft(th_white_ft)

# match-filtering the data set
xcorr2 = np.fft.irfft(sft_white*np.conj(th_white_ft))
#plt.plot(np.fft.fftshift(xcorr2))
#plt.show()



xcorr = np.fft.irfft(sft*np.fft.rfft(th[0]*window(len(th[0]),n //5)))


plt.plot(xcorr)
plt.show()

