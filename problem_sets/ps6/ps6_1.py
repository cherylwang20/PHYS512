import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy import signal
import h5py
import glob
import json

exec(open("simple_read_ligo.py").read())
f = open(r'C:\Users\wangc\PHYS512\problem_sets\ps6\LOSC_Event_tutorial\BBH_events_v3.json',)

# define the window function that is flat in the middle
# here we use the Planck-taper Window, this distribution is inspired by
# the Planck Distribution
def window(ep, N):
    x = np.arange(N)
    eN = int(ep*N)
    win = [0]*N
    win[:eN+1] = (1 + np.exp(eN/x - eN/(eN - x)))**(-1)
    win[eN: -eN] = np.ones(N - 2* eN)
    win[-eN:] = (1 + np.exp(eN/(-x + eN) - eN/x))**(-1)
    win = win[:N]
    return np.array(win)

print(len(strain_H[0]))
win = window(0.05, len(strain_H[0]))

#win = signal.tukey(len(strain_H[0]),1/5)
plt.plot(win)
plt.title('The Planck-Taper Window Function')
plt.savefig(f'Window Function.png',dpi = 300, bbox_inches = 'tight')
plt.show()


data = json.load(f)

print(json.dumps(data, indent = 4, sort_keys=True))

n  = len(strain_H[0])
nl = len(strain_L[0])
th_ft = np.abs(np.fft.rfft(th))
tl_ft = np.abs(np.fft.rfft(tl))
strain_ft = np.fft.rfft(strain_H*win)
strain_ftl = np.fft.rfft(strain_L*win)
sft = np.sum(np.abs(strain_ft),axis = 0)/4
sft_l = np.sum(np.abs(strain_ftl),axis = 0)/4


# Noise Model?
Nft_1 = np.abs(sft)**2
Nft_2 = np.abs(sft_l)**2
Nft = Nft_1.copy()
Nftl = Nft_2.copy()

# smooth out the noise model using a 1D Gaussian Filter
Nft_H = gaussian_filter1d(Nft, 1)
Nft_L = gaussian_filter1d(Nftl, 1)

plt.loglog(Nft_1, label = 'unsmoothed')
plt.loglog(Nft_H, label = 'Gaussian smoothed')
plt.legend()
plt.title('Noise Model for Hanford Detector')
plt.savefig(f'Noise_Model_H.png',dpi = 300, bbox_inches = 'tight')
plt.show()

plt.loglog(Nft_2, label = 'unsmoothed')
plt.loglog(Nft_L, label = 'Gaussian smoothed')
plt.legend()
plt.title('Noise Model for Livingston Detector')
plt.savefig(f'Noise_Model_L.png',dpi = 300, bbox_inches = 'tight')
plt.show()


sft_white = sft/np.sqrt(Nft_H)
sft_white_l = sft_l/np.sqrt(Nft_L)
th_white_ft = np.fft.rfft(th*win)/np.sqrt(Nft_H)
tl_white_ft = np.fft.rfft(tl*win)/np.sqrt(Nft_L)
th_white = np.fft.irfft(th_white_ft)
tl_white = np.fft.irfft(tl_white_ft)

# match-filtering the data set
xcorr2 = np.fft.irfft(sft_white*np.conj(th_white_ft))
xcorrl = np.fft.irfft(sft_white_l*np.conj(tl_white_ft))


# calculate the noise and SNR of each run
Noise = [0]*4
Noise_l = [0]*4
SNR = [0]*4
SNR_l = [0]*4
for i in range(4):
    Noise[i] = np.std(xcorr2[i, :-2000])
    Noise_l[i] = np.std(xcorrl[i, :-2000])
    SNR[i] = np.max(xcorr2[i]) / Noise[i]
    SNR_l[i] = np.max(xcorrl[i]) / Noise_l[i]
    print(f'#{i} GW event at Hanford has Noise of: {Noise[i]} and SNR = {SNR[i]}')
    print(f'#{i} GW event at Livingston has Noise of: {Noise_l[i]} and SNR = {SNR_l[i]}')

for i in range(4):
    plt.plot(xcorr2[i,::-1],color = 'brown')
    plt.xlim([62000,67000])
    plt.title(f'Match Filtering of #{i+1} GW event in Hanford \n with Noise = {Noise[i]:.3f} and SNR = {SNR[i]:.3f}.')
    plt.savefig(f'GW{i + 1}_H.png', dpi=300, bbox_inches='tight')
    plt.show()

for i in range(4):
    plt.plot(xcorrl[i,::-1],color = 'green')
    plt.xlim([62000,67000])
    plt.title(f'Match Filtering of #{i+1} GW event in Livingston \n with Noise = {Noise_l[i]:.3f} and SNR = {SNR_l[i]:.3f}.')
    plt.savefig(f'GW{i + 1}_L.png', dpi=300, bbox_inches='tight')
    plt.show()

xcorr = np.fft.irfft(sft*np.fft.rfft(th[0]*win))
xcorl = np.fft.irfft(sft_l*np.fft.rfft(tl[0]*win))
plt.plot(xcorr)
plt.plot(xcorl)
plt.show()

