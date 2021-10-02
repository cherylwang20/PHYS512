import time
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

matplotlib.rc('font', size = 12 , family = 'Arial')
matplotlib.rc('axes', titlesize = 12)
# define half-life in microseconds
# define converting factors from minutes, hours, days and years
mtm = 6*10**7/160
htm = 3.6*10**9/160
dtm = 8.64*10**10/160
ytm = 3.154*10**13/160

#converted to timesteps
hl = [10**21, 24*dtm, 6.7*htm, 245500*ytm, 75380*ytm, 1600*ytm,
             38235*dtm, 3.1*mtm, 26.8*mtm, 19.9*mtm, 164.3/160, 22.3*ytm,
             5.015*ytm, 138.376*dtm]


print(len(hl))

def fun(x,y,half_life = hl ):
    dydx = np.zeros(len(half_life)+1)
    dydx[0] = -y[0]/half_life[0]
    for i in range(0,len(half_life)-1):
        dydx[i+1] = y[i]/half_life[i] - y[i+1]/half_life[i+1]
    dydx[-1] = y[-2]/half_life[-2]
    return dydx

t0 = 0
t1 = 2*10.**21
y0 = [0]*(len(hl)+1)
y0[0] = 1

t_ev = np.logspace(0, np.log10(0.9*t1), num = 1000)
print(t_ev)

t_s = time.time()
# ans_rk4 = integrate.solve_ivp(fun,[t0,t1],y0,t_eval=t_ev, method = 'RK45')
ans_radau = integrate.solve_ivp(fun, [t0,t1], y0 , t_eval=t_ev, method = 'Radau')
#ans_radau = integrate.solve_ivp(fun,[t0,t1],y0 ,method = 'RK23')
t_en = time.time()

#print(ans_radau.t)
# print(ans_radau.y[0,-1],ans_radau.y[0],ans_radau.y[14],ans_radau.t)
# print(t_e - t_s)
u238 , pb206 = ans_radau.y[0],ans_radau.y[14]
th230, u234= ans_radau.y[4], ans_radau.y[3]
t = ans_radau.t
#
r1 = pb206/u238
r2 = th230/u234


plt.plot(t,r1)
plt.plot(t,u238)
plt.plot(t,pb206,'-.')
plt.savefig("upb.png",dpi=300, bbox_inches = "tight")
plt.show()

plt.plot(t[:-150],r2[:-150])
plt.plot(t,th230)
plt.plot(t, u234,'-.')
plt.yscale('log')
plt.xscale('log')
plt.savefig("thu.png",dpi=300, bbox_inches = "tight")
plt.show()