import time

import matplotlib.pyplot as plt
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

t0 = 0
t1 = 5.*10**21
y0 = [0]*(len(hl)+1)
y0[0] = 1

t_ev = np.linspace(t0,np.log(t1),400)
t_e = np.zeros(len(t_ev))
for i in range(len(t_ev)):
    t_e[i] = np.e**t_ev[i]
t_e[0] = 0

t_s = time.time()
# ans_rk4 = integrate.solve_ivp(fun,[t0,t1],y0,t_eval=t_ev, method = 'RK45')
ans_radau = integrate.solve_ivp(fun,[t0,t1],y0 ,t_eval=t_e,method = 'Radau')
t_e = time.time()

print(ans_radau.t)
# print(ans_radau.y[0,-1],ans_radau.y[0],ans_radau.y[14],ans_radau.t)
# print(t_e - t_s)
u238 , pb206 = ans_radau.y[0],ans_radau.y[14]
th230, u234= ans_radau.y[4], ans_radau.y[3]
t = ans_radau.t
#
r1 = pb206/u238
r2 = th230/u234
# print(len(t))

hfont = {'fontname':'Times'}
#
fig , (ax1, ax2) = plt.subplots(1,2 ,figsize = (12,5))
ax1.plot(t[160:-1], r1[160:-1],'r','--')
ax1.set_xlabel('Time (microseconds)',**hfont)
ax1.set_ylabel('Pb206/U238',**hfont)


ax2.plot(t[120:-40], r2[120:-40])
ax2.set_xlabel('Time (microseconds)',**hfont)
ax2.set_ylabel('Th230/U234',**hfont)

plt.savefig("Radioactive_Ratio.png", dpi = 300, bbox_inches = "tight")
plt.show()