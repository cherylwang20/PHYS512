import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower
from scipy import interpolate

dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)
step_size = np.loadtxt('step_size_i.txt')

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


# pre define a set of parameters
pars = [] # later append
#pars.append(np.asarray([60,0.02,0.1,0.05,2.00e-9,1.0]))

pars.append(np.asarray([68.74571957869041, 0.0220643823436287, 0.12071805629036651, 0.05033716294196364, 2.067806480026439e-09,
        0.9504832149645627]))
chisq = [] # later append
#model = get_spectrum([60,0.02,0.1,0.05,2.00e-9,1.0])
model = get_spectrum([68.74571957869041, 0.0220643823436287, 0.12071805629036651, 0.05033716294196364, 2.067806480026439e-09,
        0.9504832149645627])
chisq.append(get_chisq(var, model))


nstep = 1500;
step_taken = 0

def get_step(step_size):
    step = np.random.randn(6)
    stepp = step_size@step.T

    return stepp


while nstep > step_taken:
    pars_new = pars[-1] + get_step(step_size)
    print(pars[-1],chisq[-1])
    model_new = get_spectrum(pars_new)
    new_chisq = get_chisq(var,model_new)

    d_chisq = new_chisq - chisq[-1]
    prob_step = np.exp(-0.5*d_chisq)
    print(prob_step)
    accept = np.random.rand(1) < prob_step

    if accept and 0.0466< pars_new[3] < 0.0614:
        pars.append(pars_new)
        chisq.append(new_chisq)
        step_taken += 1
    print(step_taken)

print(chisq,pars[0])

output = [0]*nstep
for i in range(nstep):
    output[i] = np.append(chisq[i],pars[i])

print(output)

with open("planck_chain_tauprior.txt" , 'wb') as f:
    np.savetxt(f, output, delimiter=' ', newline='\n', header='', footer='', comments='# ')


