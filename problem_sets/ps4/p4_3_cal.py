import numpy as np
from matplotlib import pyplot as plt
dat = np.loadtxt('planck_chain_3.txt')

def DM_e(dat):
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
    return dn, uncer, err

DM_1,uncer, err = DM_e(dat)
print(f'The estimation of Dark energy is {DM_1:.3e}, with uncertainty {uncer:.3e}.')
print(err)



dat_1 = np.loadtxt('planck_chain_tauprior.txt')
dat_2 = np.loadtxt('planck_chain.txt')
dat_3 = np.loadtxt('planck_chain_2.txt')
dat_4 = np.loadtxt('planck_chain_3.txt')
dat_5 = np.loadtxt('planck_chain_tauprior_2.txt')

err_new = np.std(dat_5,axis=0)

print(err_new)
y1 = dat[:,1]
y2 = dat_2[:,1]
y3 = dat_3[:,1]
y4 = dat_4[:,1]
y5 = dat_5[:,1]

#plt.loglog(np.abs(np.fft.rfft(y2)),label = '2000 Run')
plt.loglog(np.abs(np.fft.rfft(y3)),label = '2k Run')
plt.loglog(np.abs(np.fft.rfft(y4)),color = 'green',label = '10k Run')
plt.legend()
plt.ylim(0.005,3000)
plt.savefig("MCMC_10000.png",dpi = 300, bbox_inches = 'tight')
plt.show()


plt.loglog(np.abs(np.fft.rfft(y1)),color = 'red',label = '1.5k Run')
plt.loglog(np.abs(np.fft.rfft(y5)),color = 'pink',label = '10k Run')
plt.legend()
plt.ylim(0.005,3000)
plt.savefig("MCMC_Constrain_10000.png",dpi = 300, bbox_inches = 'tight')
plt.show()

print(dat_3.shape)