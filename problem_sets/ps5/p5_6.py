import numpy as np
from matplotlib import pyplot as plt

n = 1001
x = np.arange(1,n+1)
rw = np.cumsum(np.random.randn(n))
rw_ft = np.abs(np.fft.fft(rw))
rw_ft = rw_ft/np.max(rw_ft)

#plt.plot(rw)
plt.plot(x ,rw_ft,label = 'random walk')
plt.plot(x, 1/x**2, label = '1/k Spectrum fit')
plt.xlim(-1,20)
plt.title(f'Power Spectrum of Random Walk with n = {n - 1} steps')
plt.legend()
plt.savefig(f'rw_ps.png',dpi = 300, bbox_inches = 'tight')
plt.show()

