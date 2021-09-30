import time

import numpy as np
from scipy import integrate

# define half-life in microseconds
# define converting factors from minutes, hours, days and years
mtm = 6*10**7
htm = 3.6*10**9
dtm = 8.64*10**10
ytm = 3.154*10**13

hl = [10**21, 24*dtm, 6.7*htm, 245500*ytm, 75380*ytm, 1600*ytm,
             38235*dtm, 3.1*mtm, 26.8*mtm, 19.9*mtm, 164.3, 22.3*ytm,
             5.015*ytm, 138.376*dtm]


def fun(x,y,half_life = hl ):
    dydx = np.zeros(len(half_life)+1)
    dydx[0] = -y[0]/half_life[0]
    for i in range(0,len(half_life)-1):
        dydx[i+1] = y[i]/half_life[i] - y[i+1]/half_life[i+1]
    dydx[-1] = y[-2]/half_life[-2]
    return dydx


print(len(hl))
t0 = 0
t1 = 1.7*10**21
y0 = [0]*(len(hl)+1)
y0[0] = 1

t_s = time.time()
# ans_rk4 = integrate.solve_ivp(fun,[t0,t1],y0,t_eval=t_ev, method = 'RK45')
ans_radau = integrate.solve_ivp(fun,[t0,t1],y0 , method = 'DOP853')
t_e = time.time()
print(ans_radau.y[0,-1],ans_radau.nfev)
print(t_e - t_s)