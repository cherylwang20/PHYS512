import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower
from scipy import interpolate

dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)
curv = np.loadtxt('curvature_matrix.txt').reshape(4, 2)

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];
errs = 0.5*(lowsig + highsig);

print(len(multi))

# pre define a set of parameters
pars=np.asarray([69,0.022,0.12,0.06,2.10e-9,0.95])


def get_spectrum(pars,lmax=3000):
    #print('pars are ',pars)
    H0=pars[0]
    ombh2=pars[1]
    omch2=pars[2]
    tau=pars[3]
    As=pars[4]
    ns=pars[5]
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]
    tt = tt[2:]
    return tt[:len(var)]

def deriv(f, x, n):
    dx = 0.01*x
    x1 = x + 2 * dx
    x2 = x + dx
    x3 = x - dx
    x4 = x - 2 * dx
    return (f(x4,n) - 8 * f(x3,n) + 8 * f(x2,n) - f(x1,n)) / (12 * dx)

def ff(x,n):
    pars = np.asarray([60,0.02,0.1,0.05,2.00e-9,1.0])
    pars[n] = x
    print(pars)
    v = get_spectrum(pars)
    return v

der = np.zeros([len(multi), len(pars)])

#create the dervative matrix
for n in range(len(pars)):
    der[:,n] = deriv(ff, pars[n],n)

def step_size():

    return