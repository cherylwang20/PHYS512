import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

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
def ratfit(y,x,n,m):
    npt=len(x)
    assert(len(y)==npt)
    assert(n>=0)
    assert(m>=0)
    assert(n+ 1+m==npt)

    top_mat=np.empty([npt,n+1])
    bot_mat=np.empty([npt,m])
    for i in range(n+1):
        top_mat[:,i]=x**i
    for i in range(m):
        bot_mat[:,i]=y*x**(i+1)
    mat=np.hstack([top_mat,-bot_mat])
    print(mat)
    pars=np.linalg.inv(mat)@y
    p=pars[:n+1]
    q=pars[n+1:]
    return mat,p,q

def rateval(x,p,q):
    top=0
    for i,par in enumerate(p):
        top=top+par*x**i
    bot=1
    for i,par in enumerate(q):
        bot=bot+par*x**(i+1)
    return top/bot

#n = 9
#m = 10 # why is it that the larger the value the greater the error?
xr = np.arange(-10,10)
yr = np.cos(xr)
m = len(yr) //2
n = len(yr) - m - 1
mat, p, q = ratfit(yr, xr, n, m)
yy_r = rateval(xx, p, q)

plt.plot(xx,yy_r)
#plt.plot(xx,yy_p,'red')
#plt.plot(xx,yy_c,'-.')
plt.show()

y_cos = np.cos(xx)
print('error through polynomial interpolation of cosine function is', np.std(yy_p-y_cos))
print('error through Cubic Spline interpolation of cosine function is', np.std(yy_c-y_cos))
print('error through Rational Function interpolation of cosine function is', np.std(yy_r-y_cos))

# for the Lorentz function
# complete the three methods for this function

# the rational function fit
# p_l = [1,2]
# q_l = [0, 1]
# x_l = np.linspace(-1, 1, 1001)
# y_l = rat_eval(p_l,q_l,x_l)
#
# y_lr = 1/(1+x_l**2)
#
# err = np.std(y_lr - y_l)
# print(err)


