import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import json

exec(open("simple_read_ligo.py").read())
f = open(r'C:\Users\wangc\PHYS512\problem_sets\ps6\LOSC_Event_tutorial\BBH_events_v3.json',)

# define the window function that is flat in the middle
def window(n):
    x = np.linspace(-np.pi, np.pi, n)
    win = 0.5 - 0.5*np.cos(2*np.pi*x)
    return win

data = json.load(f)

th_ft = np.abs(np.fft.rfft(th))
tl_ft = np.abs(np.fft.rfft(tl))
strain_ft = np.abs(np.fft.rfft(strain*window(len(strain))))




plt.plot(th)
#plt.plot(tl)
plt.show()


plt.loglog(strain_ft)
plt.show()