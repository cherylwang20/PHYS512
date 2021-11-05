import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import quad
import scipy
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
    win[:eN] = (1 + np.exp(eN/x - eN/(eN - x)))**(-1)
    win[eN: -eN] = np.ones(N - 2* eN)
    win[-eN:] = (1 + np.exp(eN/(-x + eN) - eN/x))**(-1)
    win = win[:N]
    return np.array(win)

print(len(strain_H[0]))
win = window(0.1, len(strain_H[0]))

#win = signal.tukey(len(strain_H[0]),1/5)
plt.plot(win)
plt.title('The Planck-Taper Window Function')
plt.savefig(f'Window Function.png',dpi = 300, bbox_inches = 'tight')
plt.show()


data = json.load(f)

print(json.dumps(data, indent = 4, sort_keys=True))

n  = len(strain_H[0])
nl = len(strain_L[0])
strain_ft = np.fft.rfft(strain_H*win)
strain_ftl = np.fft.rfft(strain_L*win)


# Noise Model?
Nft_1 = np.abs(strain_ft)**2
Nft_2 = np.abs(strain_ftl)**2
Nft = Nft_1.copy()
Nftl = Nft_2.copy()

Nft_H = [0]*4
Nft_L = [0]*4
# smooth out the noise model using a 1D Gaussian Filter
#Nft_H = gaussian_filter1d(Nft, 1)
for i in range(4):
    Nft_H[i] = scipy.signal.medfilt(Nft[i], 11)
#Nft_L = gaussian_filter1d(Nftl, 1)
    Nft_L[i] = scipy.signal.medfilt(Nftl[i], 11)

#Nft_L = gaussian_filter1d(Nftl, 1)
# for i in range(4):
#     plt.loglog(Nft_1[i], label = f'unsmoothed {i + 1}')
#     plt.loglog(Nft_H[i], label = f'Gaussian smoothed {i + 1}')
#     plt.legend()
#     plt.title(f'Noise Model for Hanford Detector {i + 1}')
#     plt.savefig(f'Noise_Model_H{i+1}.png',dpi = 300, bbox_inches = 'tight')
#     plt.show()
#
#     plt.loglog(Nft_2[i], label = 'unsmoothed')
#     plt.loglog(Nft_L[i], label = 'Gaussian smoothed')
#     plt.legend()
#     plt.title(f'Noise Model for Livingston Detector {i + 1}')
#     plt.savefig(f'Noise_Model_L{i+1}.png',dpi = 300, bbox_inches = 'tight')
#     plt.show()
#whi

# whitening the function
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    freqs1 = np.linspace(0,2048.,Nt/2+1)

    # whitening: transform to freq domain, divide by asd, then transform back,
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    norm = 1./np.sqrt(1./(dt*2))
    white_hf = hf / np.sqrt(interp_psd(freqs)) * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

#sft_white = whiten(strain_ft, scipy.integrate.interp1d() )
sft_white = strain_ft /np.sqrt(Nft_H)
sft_white_l = strain_ftl/np.sqrt(Nft_L)
th_white_ft = np.fft.rfft(th*win)/np.sqrt(Nft_H)
tl_white_ft = np.fft.rfft(tl*win)/np.sqrt(Nft_L)
th_white = np.fft.irfft(th_white_ft)
tl_white = np.fft.irfft(tl_white_ft)

# match-filtering the data set
xcorr2 = np.fft.irfft(sft_white*np.conj(th_white_ft))
xcorrl = np.fft.irfft(sft_white_l*np.conj(tl_white_ft))

# calculate the frequency to achieve half of the weight
fxch = np.abs(np.fft.rfft(xcorr2))
ftot = np.sum(fxch, axis = 1)
for j in range(len(ftot)):
    for i in range(len(fxch[j])):
        fhalf = np.sum(fxch[j,:i])
        if fhalf > ftot[j]/2:
            print(f'The half weight frequency for the #{j + 1} Event at Hanford is {i:.3e} Hz')
            break

fxcl = np.abs(np.fft.rfft(xcorrl))
ftotl = np.sum(fxcl, axis = 1)
for j in range(len(ftotl)):
    for i in range(len(fxcl[j])):
        fhalfl = np.sum(fxcl[j,:i])
        if fhalfl > ftotl[j]/2:
            print(f'The half weight frequency for the #{j + 1} Event at Livingston is {i:.3e} Hz')
            break
plt.loglog(fxch[0])
plt.show()


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
    print(f'#{i + 1} GW event at Hanford has Noise of: {Noise[i]} and SNR = {SNR[i]}')
    print(f'#{i + 1} GW event at Livingston has Noise of: {Noise_l[i]} and SNR = {SNR_l[i]}')

for i in range(4):
    plt.plot(np.fft.fftshift(xcorr2[i]),color = 'brown')
    #plt.xlim([62000,67000])
    plt.title(f'Match Filtering of #{i+1} GW event in Hanford \n with Noise = {Noise[i]:.3f} and SNR = {SNR[i]:.3f}.')
    plt.savefig(f'GW{i + 1}_H.png', dpi=300, bbox_inches='tight')
    plt.show()

for i in range(4):
    plt.plot(np.fft.fftshift(xcorrl[i]),color = 'green')
    #plt.xlim([62000,67000])
    plt.title(f'Match Filtering of #{i+1} GW event in Livingston \n with Noise = {Noise_l[i]:.3f} and SNR = {SNR_l[i]:.3f}.')
    plt.savefig(f'GW{i + 1}_L.png', dpi=300, bbox_inches='tight')
    plt.show()


sigma_h = np.sqrt(np.abs(np.fft.irfft(th_white_ft*np.conj(th_white_ft))))
sigma_l = np.sqrt(np.abs(np.fft.irfft(tl_white_ft*np.conj(tl_white_ft))))

Noise_ah = [0]*4
Noise_al = [0]*4
SNR_ah= [0]*4
SNR_al = [0]*4

for i in range(4):
    Noise_ah[i] = np.std(sigma_h[i, :-2000])
    Noise_al[i] = np.std(sigma_l[i, :-2000])
    SNR_ah[i] = np.max(sigma_h[i]) / Noise_ah[i]
    SNR_al[i] = np.max(sigma_l[i]) / Noise_al[i]
    print(f'#{i + 1} Noise Model at Hanford has Noise of: {Noise_ah[i]} and SNR = {SNR_ah[i]}')
    print(f'#{i + 1} Noise Model at Livingston has Noise of: {Noise_al[i]} and SNR = {SNR_al[i]}')

for i in range(4):
    plt.plot(np.fft.fftshift(sigma_h[1]))
    plt.title(f'Analytic Noise Model for Hanford #{i + 1}')
    plt.show()

    plt.plot(np.fft.fftshift(sigma_l[1]))
    plt.title(f'Analytic Noise Model for Livingston #{i + 1}')
    plt.show()



# two LIGO detectors are 3000 km apart
print(utc_H)
print(utc_L)