import numpy as np
from matplotlib import pyplot as plt

n = 10001
x = np.arange(1,n+1)
y = 1/x**2
rw = np.cumsum(np.random.randn(n))
rw_ft = np.abs(np.fft.fft(rw))
rw_ft = rw_ft/np.max(rw_ft)

win =  0.5 - 0.5*np.cos(2*np.pi*x/n)
rw_win = rw*win
rw_win_ft = np.abs(np.fft.fft(rw_win))
rw_win_ft = rw_win_ft/np.max(rw_win_ft)




#plt.plot(rw)
plt.plot(x ,rw_ft,label = 'random walk w/o window')
plt.plot(x, rw_win_ft, label = 'random walk w/ window')
plt.plot(x, y, label = '1/k Spectrum fit')
plt.xlim(-1,20)
plt.title(f'Power Spectrum of Random Walk with n = {n - 1} steps')
plt.legend()
#plt.savefig(f'rw_ps.png',dpi = 300, bbox_inches = 'tight')
plt.show()

chisq_win = np.sum((rw_win_ft[:21] - y[:21])**2/y[:21])
chisq = np.sum((rw_ft[:21] - y[:21])**2/y[:21])
print(f'A 1/x^2 model gives a chisq of {chisq} for unwindowed transform, which is desirable.')
print(f'A 1/x^2 model gives a chisq of {chisq_win} for windowed transform, which is desirable.')


tt = -3*x
tf = np.abs(np.fft.fft(tt))/np.sum(np.abs(np.fft.fft(tt)))
plt.plot(x,tf)
plt.xlim(-1,20)
plt.show()