import numpy as np
from scipy import interpolate

# polynomial interpolation
x1 = np.linspace(-np.pi,np.pi, 1001)
dx = x1[1] - x1[0]
y1 = np.cos(x1)

xx = np.linspace(-np.pi/2, np.pi/2, 1001)
yy_p = np.empty(len(xx))
for i in range(len(xx)):
    ind = (xx[i] - x1[0])/dx
    ind = int(np.floor(ind))
    x_use = x1[ind-1:ind+3]
    y_use = y1[ind-1:ind+3]
    p = np.polyfit(x_use,y_use,3)
    yy_p[i] = np.polyval(p,xx[i])

# cubic spline fit
spln = interpolate.splrep(x1,y1)
yy_c = interpolate.splev(xx,spln)

# rational function fit
def rat_fit(x,y,n,m):
    assert (len(x) == n + m - 1)
    assert (len(y) == len(x))
    mat = np.zeros([n+m - 1, n+m - 1])
    for i in range(n):
        mat[:,i] = x**i
    for i in range(1, m):
        mat[:, i - 1+n] = -y*x**i
    pars = np.dot(np.linalg.inv(mat),y)
    p = pars[:n]
    q = pars[n:]
    return p,q

def rat_eval(p,q,x):
    top = 0
    for i in range(len(p)):
        top = top + p[i]*x**i
    bot = 1
    for i in range(len(q)):
        bot = bot + q[i]*x**(i+1)
    return top/bot

n = 4
m = 5 # why is it that the larger the value the greater the error?
xr = np.linspace(-np.pi/2, np.pi/2, n+m-1)
yr = np.cos(xr)
p, q = rat_fit(xr, yr, n, m)
yy_r = rat_eval(p, q, xx)


y_cos = np.cos(xx)
print('error through polynomial interpolation of cosine function is', np.std(yy_p-y_cos))
print('error through Cubic Spline interpolation of cosine function is', np.std(yy_c-y_cos))
print('error through Rational Function interpolation of cosine function is', np.std(yy_r-y_cos))

# for the Lorentz function
# complete the three methods for this function

# the rational function fit
p_l = [1,2]
q_l = [0, 1]
x_l = np.linspace(-1, 1, 1001)
y_l = rat_eval(p_l,q_l,x_l)

y_lr = 1/(1+x_l**2)

err = np.std(y_lr - y_l)
print(err)


