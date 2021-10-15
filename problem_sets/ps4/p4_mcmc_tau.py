import numpy as np
import camb
from camb import model, initialpower
from scipy import interpolate

data = np.loadtxt('planck_chain.txt')
tau = 0.0540
sigma_tau = 2.*0.0074

H0 = data[:,1]
bd = data[:,2]
dmd = data[:,3]
ti = data[:,4]
As = data[:,5]
ns = data[:,6]

p_i = np.exp(-0.5 * ((ti - tau)/sigma_tau)**2)
nor = np.sum(p_i)

W = [0]*6
for i in range(6):
    W[i] = np.sum(data[:,i+1]*p_i)/nor

curv = np.loadtxt('curvature_matrix.txt')
step_size = np.linalg.cholesky(np.linalg.inv(curv))

pars = W
print(W)

def get_spectrum(pars,lmax=3000):
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
    tt=cmb[:,0]    #you could return the full power spectrum here if you wanted to do say EE
    tt = tt[2:]
    return tt[:len(p_i)]

# finding the derivative matrix
def deriv(f, x, n):
    dx = 0.01*x
    x1 = x + 2 * dx
    x2 = x + dx
    x3 = x - dx
    x4 = x - 2 * dx
    return (f(x4,n) - 8 * f(x3,n) + 8 * f(x2,n) - f(x1,n)) / (12 * dx)

def ff(x,n):
    pars = [68.74571957869041, 0.0220643823436287, 0.12071805629036651, 0.05033716294196364, 2.067806480026439e-09, 0.9504832149645627]
    pars[n] = x
    print(pars)
    v = get_spectrum(pars)
    return v

der = np.zeros([len(H0), len(pars)])

#create the dervative matrix
for n in range(len(pars)):
    der[:,n] = deriv(ff, pars[n],n)

lhs=der.T@der
step_size = np.linalg.cholesky(np.linalg.inv(lhs))


with open("step_size_i.txt" , 'wb') as f:
    np.savetxt(f, step_size, delimiter=' ', newline='\n', header='', footer='', comments='# ')
