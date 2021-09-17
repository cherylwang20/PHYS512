import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# polynomial interpolation
# we take in a function with a initial range of data points
# we also take in the final points/interval that we want would like to interpolate
def poly(fun,x1,dx,xx):
    # create an empty array for final interpolation
    yy_p = np.empty(len(xx))
    # evaluate the values of the function at x1s
    y1 = fun(x1)
    # repeat the process for the desire array xx and append it to yy_p
    for i in range(len(xx)):
        ind = (xx[i] - x1[0])/dx
        ind = int(np.floor(ind))
        x_use = x1[ind-1:ind+3]
        y_use = y1[ind-1:ind+3]
        p = np.polyfit(x_use,y_use,3)
        yy_p[i] = np.polyval(p,xx[i])
        err = np.std(yy_p - fun(xx)) #compare the interpolated data with the actual function
    return yy_p,err

# cubic spline fit
# we intake the function, the initial data points/range
# also take in the final desired interpolated data
def cubint(fun,x1,xx):
    y1 = fun(x1)
    spln = interpolate.splrep(x1,y1)
    yy_c = interpolate.splev(xx,spln)
    err = np.std(yy_c - fun(xx)) #evalute the interpolated value with the actual function
    return yy_c,err

# rational function fit
def ratfit(y,x,n,m, full = False):
    # intake the length of x and assert out matrix size is correct
    npt=len(x)
    assert(len(y)==npt)
    assert(n>=0)
    assert(m>=0)
    assert(n+ 1+m==npt)

    # create an empty matrix
    top_mat=np.empty([npt,n+1])
    bot_mat=np.empty([npt,m])
    # append the matrix with the desired value
    for i in range(n+1):
        top_mat[:,i]=x**i
    for i in range(m):
        bot_mat[:,i]=y*x**(i+1)
    # merge the top and bottom row
    mat=np.hstack([top_mat,-bot_mat])

    # create a choice for user to use either pinv and inv
    if full == True:
        pars =np.linalg.pinv(mat)@y
    else:
        pars=np.linalg.inv(mat)@y

    # generate the coefficient and p and q
    p=pars[:n+1]
    q=pars[n+1:]
    return mat,p,q

def rateval(x,p,q):
    # using the cofficients generated form ratfit, compute the top and bottom
    # polynomial functions for rational function fit
    print(p,q)
    top=0
    for i,par in enumerate(p):
        top=top+par*x**i
    bot=1
    for i,par in enumerate(q):
        bot=bot+par*x**(i+1)
    return top/bot

# intake the function, the initial xr values and the final desired value
# for interpolation
def ratint(fun,xr,xx,full = False):
    # create data points of the function
    yr = fun(xr)
    # define our m and n of the matrix
    m = len(yr) //2
    n = len(yr) - m - 1
    # the user may use eithr pinv or inv fit
    if full == True:
        mat, p, q = ratfit(yr, xr, n, m, True) # using pinv fit
    else:
        mat, p, q = ratfit(yr, xr, n, m) # using inv fit
    yy_r = rateval(xx, p, q)
    err = np.std(yy_r - fun(xx)) # evaluate the interpolated data with actual function
    return yy_r,err

# define the lorentz function
def lor(x):
    return 1/(1+x**2)

if __name__ == '__main__':
    # the first function to evaluate in cosine,
    x1 = np.linspace(-np.pi, np.pi, 1001) # choose the initial data set to interpolate
    dx = x1[1] - x1[0] # define the desired interval
    xx = np.linspace(-np.pi / 2, np.pi / 2, 1001) # choose the final desired data out of interpolation

    polycos = poly(np.cos, x1, dx, xx)
    print('error through polynomial interpolation of,', np.cos, ' is', f" {polycos[1]: .3e}")

    cubcos = cubint(np.cos, x1, xx)
    print('error through Cubic Spline interpolation of ,', np.cos, ' is', f" {cubcos[1]: .3e}")


    xr = np.arange(-10,10)
    ratcos = ratint(np.cos,xr,xx)
    print('error through Rational Function interpolation of  ,',np.cos,' is', f" {ratcos[1]: .3e}" )



    # for the Lorentz function, same as the cosine function
    xl = np.linspace(-5,5,1001)
    dxl = xl[1] - xl[0]
    xxl = np.linspace(-1,1,1001)
    polylor = poly(lor,xl,dxl,xxl)
    print('error through polynomial interpolation of the lorenz function is', f" {polylor[1]: .3e}")

    cublor = cubint(lor,xl,xxl)
    print('error through Cubic Spline interpolation of the lorentz function is', f" {cublor[1]: .3e}")

    xrl = np.arange(-3,2)
    ratlor = ratint(lor,xrl,xxl)
    print('error through Rational Function interpolation the lorenz function is', f" {ratlor[1]: .3e}")

    xrl_2 = np.arange(-5,5)
    ratlor_2 = ratint(lor, xrl_2, xxl)
    print('error through Rational Function interpolation the lorenz function with n = 4, m = 5 is', f" {ratlor_2[1]: .3e}")

    ratlor_2p = ratint(lor, xrl_2, xxl, full = True)
    print('if we use linalg.pinv instead of linalg.inv, we have the error fit for Lorenz function through Rational function interpolation as', f" {ratlor_2p[1]: .3e}" )
