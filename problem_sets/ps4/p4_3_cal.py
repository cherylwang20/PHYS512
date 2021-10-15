import numpy as np

dat = np.loadtxt('planck_chain.txt')

#assuming errors are gaussian
err = np.std(dat,axis=0)
mean = np.mean(dat,axis=0)

H0 = mean[1]
h = H0/100.
bd = mean[2]
dmd = mean[3]

dn = 1 - bd/h**2 - dmd/h**2
sighsq = err[1]*2/100.
sigbd = err[2]
sigdmd = err[3]

#compute error using propagation of uncertainty
uncer1 = bd/h**2*np.sqrt((sigbd/bd)**2 + (sighsq/h**2)**2 - 2*(sigbd*sighsq)/(bd*h**2))
uncer2 = dmd/h**2*np.sqrt((sigdmd/dmd)**2 + (sighsq/h**2)**2 - 2*(sigdmd*sighsq)/(dmd*h**2))
uncer =uncer2 + uncer1

print(err)
print(f'The estimation of Dark energy is {dn:.3e}, with uncertainty {uncer:.3e}.')
