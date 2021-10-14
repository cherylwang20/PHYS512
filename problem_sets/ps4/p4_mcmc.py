import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower
from scipy import interpolate

dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)
curv = np.loadtxt('curvature_matrix.txt')
der = np.loadtxt('deriv_matrix.txt')

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];
errs = 0.5*(lowsig + highsig);



def get_chisq(data,model):
    chisq = np.sum((data-model)**2/errs**2)
    return chisq

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

# pre define a set of parameters
pars = [] # later append
#pars.append(np.asarray([60,0.02,0.1,0.05,2.00e-9,1.0]))
pars.append(np.asarray([65,0.022, 0.12, 0.06, 2.10e-9,0.95 ]))
chisq = [] # later append
#model = get_spectrum([60,0.02,0.1,0.05,2.00e-9,1.0])
model = get_spectrum([65,0.02, 0.1, 0.05, 2.00e-9,1 ])
chisq.append(get_chisq(var, model))
print(get_chisq(var, model))


nstep = 10;
step_size = np.linalg.inv(curv)@der.T@((var - model)/3)
print(step_size)
step_taken = 0

def get_step(step_size):
    step = np.random.randn(len(step_size))*step_size
    return step

while nstep > step_taken:
    pars_new = pars[-1] + get_step(step_size)
    print(pars[-1],chisq[-1])
    model_new = get_spectrum(pars_new)
    new_chisq = get_chisq(var,model_new)

    d_chisq = new_chisq - chisq[-1]
    prob_step = np.exp(-0.5*d_chisq)
    print(prob_step)
    accept = np.random.rand(1) < prob_step

    if accept:
        pars.append(pars_new)
        chisq.append(new_chisq)
        step_taken += 1
    print(step_taken)
    par_errs = np.sqrt(np.diag(np.linalg.inv(curv)))



print(pars)

with open("planck_chain.txt" , 'wb') as f:
    np.savetxt(f, pars, delimiter=' ', newline='\n', header='', footer='', comments='# ')


