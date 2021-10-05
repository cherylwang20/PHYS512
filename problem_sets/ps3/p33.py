import numpy as np

dat = np.loadtxt('dish_zenith.txt')

# position coordinate in mm
x = dat[:,0]
y = dat[:,1]
z = dat[:,2]

A = np.zeros([len(x),4])

A[:,0] =  x**2 + y**2
A[:,1] = -x
A[:,2] = -y
A[:,3] = [1.]*len(x)

A = np.matrix(A)

d = np.matrix(z).transpose()

u,s,v = np.linalg.svd(A,False)
sinv = np.matrix(np.diag(1.0/s))
fitp = v.transpose()*sinv*(u.transpose()*d)

fitp = np.squeeze(np.asarray(fitp))
print(fitp)

q, b, c, p = fitp
a = q
x0 = b/2/a
y0 = c/2/a
z0 = p - a*x0**2 - a*y0**2
print(a,x0,y0,z0)

zn = a*((x-x0)**2 + (y - y0)**2) + z0

noise = np.std(abs(z - zn))
N = np.eye(len(x))*noise**2
Ninv = np.eye(len(x))*noise**(-2)
mat=A.T@Ninv@A
errs = np.linalg.inv(mat)
err = np.sqrt(np.diag(errs))

derrs = A@errs@A.T
model_sig = np.sqrt(np.diag(derrs))

print(noise)
print('The uncertainty of a is', err[0])

