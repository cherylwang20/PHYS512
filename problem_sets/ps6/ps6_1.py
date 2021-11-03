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
th_ft = np.abs(np.fft.rfft(th))
tl_ft = np.abs(np.fft.rfft(tl))
strain_ft = np.fft.rfft(strain_H*win)
sft = np.sum(np.abs(strain_ft),axis = 0)/4

print(sft)

# for i in range(4):
#     plt.loglog(np.abs(strain_ft[i]))
# plt.show()

# Noise Model?
Nft_1 = np.abs(sft)**2
Nft = Nft_1.copy()
# smooth out the noise model:
for i in range(10):
    Nft_H = (Nft + np.roll(Nft, 1) + np.roll(Nft, -1))/3

# plt.loglog(Nft_1, label = 'unsmoothed')
# plt.loglog(Nft_H, label = 'smoothed')
# plt.legend()
# plt.savefig(f'Noise_Model_H.png',dpi = 300, bbox_inches = 'tight')
# plt.show()


sft_white = sft/np.sqrt(Nft_H)
th_white_ft = np.fft.rfft(th*win)/np.sqrt(Nft_H)
th_white = np.fft.irfft(th_white_ft)

# match-filtering the data set
xcorr2 = np.fft.irfft(sft_white*np.conj(th_white_ft))
#plt.plot(np.fft.fftshift(xcorr2))

Noise = [0]*4
SNR = [0]*4
for i in range(4):
    Noise[i] = np.std(xcorr2[i, :-2000])
    SNR[i] = np.max(xcorr2[i]) / Noise[i]
    print(f'#{i} GW event has Noise of: {Noise[i]} and SNR = {SNR[i]}')

for i in range(4):
    plt.plot(xcorr2[i,::-1],color = 'brown')
    plt.xlim([56000,67000])
    plt.title(f'Match Filtering of #{i+1} GW event in Hanford \n with Noise = {Noise[i]:.3f} and SNR = {SNR[i]:.3f}.')
    plt.savefig(f'GW{i + 1}_H.png', dpi=300, bbox_inches='tight')
    plt.show()

# Noise = np.std(xcorr2[0,:-2000])
# print(np.std(xcorr2[0,2000:]),np.std(xcorr2[0,:-2000]))
# SNR = np.max(xcorr2[0])/Noise




xcorr = np.fft.irfft(sft*np.fft.rfft(th[0]*window(len(th[0]),n //5)))
plt.plot(xcorr)
plt.show()

