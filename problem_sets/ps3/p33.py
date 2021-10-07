import matplotlib.pyplot as plt
import numpy as np

dat = np.loadtxt('dish_zenith.txt')

# extract position coordinate in mm
x = dat[:,0]
y = dat[:,1]
z = dat[:,2]

# create a matrix and define its parameters
A = np.zeros([len(x),4])
A[:,0] =  x**2 + y**2
A[:,1] = -x
A[:,2] = -y
A[:,3] = [1.]*len(x)
A = np.matrix(A)

d = np.matrix(z).transpose()

# carry out the fit using SVD
u,s,v = np.linalg.svd(A,False)
sinv = np.matrix(np.diag(1.0/s))
fitp = v.transpose()*sinv*(u.transpose()*d)

fitp = np.squeeze(np.asarray(fitp))
print(f'The best fit parameters are: p = {fitp[0]:.3e}, b = {fitp[1]:.3e}, c = {fitp[2]:.3e}, d = {fitp[3]:.3e}')

# convert the matrix element into the parameters of desire
q, b, c, p = fitp
a = q
x0 = b/2/a
y0 = c/2/a
z0 = p - a*x0**2 - a*y0**2

print(f'The best fit parameters are: a = {a:.3e}, x_0 = {x0:.3f}, y_0 = {y0:.3f}, z_0 = {z0:.3f}')


# compute our z
zn = a*((x-x0)**2 + (y - y0)**2) + z0

# compare with the original value of z and acquire the noise
noise = np.std(z - zn)

# plot the actual noise of the system
plt.plot((z - zn))
plt.ylabel('Actual Noise')
plt.savefig("Noise.png",dpi=300, bbox_inches = "tight")
plt.show()

# compute the error bar on the data set we have
N = np.eye(len(x))*noise**2
Ninv = np.eye(len(x))*noise**(-2)
mat=A.T@Ninv@A
errs = np.linalg.inv(mat)
err = np.sqrt(np.diag(errs)) # error estimation on each parameter

derrs = A@errs@A.T
model_sig = np.sqrt(np.diag(derrs)) # the overall noise in the data

print(f'The noise of the system is: {noise:.3f}. The overall model sigma is: {np.std(model_sig):.3f}.')
print(f'The uncertainty of a is: {err[0]:.3e}')

# calculating the focal length, derivation in doc
r = np.sqrt(-z0/a)+y0
f = abs(r**2/4/z0)/1000 # in m

# error bar
err_bar =  (r/2/z0*(noise**2) - r**2/4/z0**2*(noise**2))/1000

print(f'The focal length is:  {f:.4f} +/- {-err_bar:.3f} m')



