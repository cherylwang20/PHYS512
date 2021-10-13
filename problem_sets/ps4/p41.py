import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower
from scipy import interpolate


dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];
errs = 0.5*(lowsig + highsig);

print(len(multi))

# pre define a set of parameters
pars=np.asarray([60,0.02,0.1,0.05,2.00e-9,1.0])


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
    return tt[:len(multi)]

def ff(x,n):
    pars[n] = x
    v = get_spectrum(pars)
    return v

# numerial derivaives
def ndiff(fun,x,n):
    dx = d(fun,x,n)
    print(dx)
    deriv = (fun(x + dx,n) - fun(x - dx,n))/(2*dx)
    return deriv

def d(f, x,n):
    eps = 7 * 10 ** (-18)
    # choose the optimal h with trials
    h = eps ** (1 / 2)
    # derive the function for the third derivative use to estimated dx
    func = (f(x+2*h,n) - 3*f(x+h,n) + 3*f(x,n) - f(x- h,n))/ h ** 3
    # evaluate the optimal delta value for each input x and function f
    delt = (3 * eps * abs(f(x,n)) / abs(func)) ** (1 / 3)
    return np.std(delt)

derivs = np.zeros([len(multi), len(pars)])

for n in range(len(pars)):
    derivs[:,n] = ndiff(ff,pars[n],n)

print(derivs)