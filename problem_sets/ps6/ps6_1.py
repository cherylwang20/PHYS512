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

n  = len(strain)
win = window(len(strain), n//5)
th_ft = np.abs(np.fft.rfft(th))
tl_ft = np.abs(np.fft.rfft(tl))
strain_ft = np.fft.rfft(strain*win)

# Noise Model?
Nft = np.abs(strain_ft)**2
# smooth out the noise model:
for i in range(10):
    Nft = (Nft + np.roll(Nft, 1) + np.roll(Nft, -1))/3
sft_white = strain_ft/np.sqrt(Nft)
th_white_ft = np.fft.rfft(th*win)/np.sqrt(Nft)
th_white = np.fft.irfft(th_white_ft)

# match-filtering the data set
xcorr2 = np.fft.irfft(sft_white*np.conj(th_white_ft))
plt.plot(np.fft.fftshift(xcorr2))
plt.show()



plt.loglog(Nft)
plt.show()

xcorr = np.fft.irfft(strain_ft*np.fft.rfft(th*window(len(th),n //5)))


plt.plot(xcorr)
#plt.plot(tl)
plt.show()


plt.loglog(np.abs(strain_ft))
plt.show()