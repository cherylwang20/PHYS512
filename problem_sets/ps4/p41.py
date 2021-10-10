import numpy as np

dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];
