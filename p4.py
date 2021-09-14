import numpy as np
from scipy import interpolate

x1 = np.linspace(-np.pi,np.pi, 1001)
dx = x1[1] - x1[0]
y1 = np.cos(x1)

# polynomial interpolation
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



y_cos = np.cos(xx)
print('error through polynomial interpolation of cosine function is', np.std(yy_p-y_cos))
print('error through Cubic Spline interpolation of cosine function is', np.std(yy_c-y_cos))

