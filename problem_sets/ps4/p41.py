import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower


dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];

