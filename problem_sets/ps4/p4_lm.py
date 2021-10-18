import numpy as np
from matplotlib import pyplot as plt
import camb
from camb import model, initialpower
from scipy import interpolate


dat = np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt',skiprows=1)

multi = dat[:,0]; var = dat[:,1];
lowsig = dat[:,2]; highsig = dat[:,3];
errs = 0.5*(lowsig + highsig);

pars=np.asarray([65,0.021,0.11,0.055,2.05e-9,0.97])
#pars_new = np.asarray([69,0.022, 0.12, 0.06, 2.10e-9,0.95 ])

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
    return tt[:len(var)]

# finding the derivative matrix
def deriv(f, x, n):
    dx = 0.01*x
    x1 = x + 2 * dx
    x2 = x + dx
    x3 = x - dx
    x4 = x - 2 * dx
    return (f(x4,n) - 8 * f(x3,n) + 8 * f(x2,n) - f(x1,n)) / (12 * dx)

def ff(x,n):
    pars = [65,0.021,0.11,0.055,2.05e-9,0.97]
    pars[n] = x
    print(pars)
    v = get_spectrum(pars)
    return v

der = np.zeros([len(multi), len(pars)])

#create the dervative matrix
for n in range(len(pars)):
    der[:,n] = deriv(ff, pars[n],n)

with open("deriv_matrix.txt" , 'wb') as f:
    np.savetxt(f, der, delimiter=' ', newline='\n', header='', footer='', comments='# ')

def get_matrices(pars, fun, der, y):
    model = fun(pars)
    derivs = der
    r = y - model
    lhs = derivs.T @ derivs
    rhs = derivs.T @ r
    chisq = np.sum((r/errs) ** 2)
    return chisq, lhs, rhs

def linv(mat, lamda):
    mat = mat + lamda * np.diag(np.diag(mat))
    return np.linalg.inv(mat)


def fit_lm_clean(m, fun, der, y, niter=20, chitol=0.01):
    lamda = 0
    chisq, lhs, rhs = get_matrices(m, fun, der, y)
    for i in range(niter):
        lhs_inv = linv(lhs, lamda)
        dm = lhs_inv @ rhs
        chisq_new, lhs_new, rhs_new = get_matrices(m + dm, fun, der, y)
        if chisq_new < chisq:
            if lamda == 0:
                if (np.abs(chisq - chisq_new) < chitol):
                    print(np.abs(chisq - chisq_new))
                    print('Converged after ', i, ' iterations of LM')
                    return m + dm
            chisq = chisq_new
            lhs = lhs_new
            rhs = rhs_new
            m = m + dm
            lamda = update_lamda(lamda, True)

        else:
            lamda = update_lamda(lamda, False)
        print('on iteration ', i, ' chisq is ', chisq, ' with step ', dm, ' and lamda ', lamda)
        print(m)
    par_errs = np.sqrt(np.diag(np.linalg.inv(lhs)))
    return m, par_errs, lhs

def update_lamda(lamda,success):
    if success:
        lamda=lamda/1.5
        if lamda<0.5:
            lamda=0
    else:
        if lamda==0:
            lamda=1
        else:
            lamda=lamda*1.5**2
    return lamda

fit = fit_lm_clean(pars,get_spectrum,der,var)
#fit_new = fit_lm_clean(pars_new,get_spectrum,der,var)

curv = fit[2]
with open("curvature_matrix.txt" , 'wb') as f:
    np.savetxt(f, curv, delimiter=' ', newline='\n', header='', footer='', comments='# ')


output = [fit[0],fit[1]]
with open("planck_fit_params_lm.txt" , 'wb') as f:
    np.savetxt(f, output, delimiter=' ', newline='\n', header='', footer='', comments='# ')


# output_new = [fit_new[0],fit_new[1]]
# with open("planck_fit_params_lm_new.txt" , 'wb') as f:
#     np.savetxt(f, output_new, delimiter=' ', newline='\n', header='', footer='', comments='# ')
