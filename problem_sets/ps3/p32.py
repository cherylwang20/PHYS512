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
             3.8235*dtm, 3.1*mtm, 26.8*mtm, 19.9*mtm, 164.3/160, 22.3*ytm,
             5.015*ytm, 138.376*dtm]


print(len(hl))

# define the function that runs through the lifetime of each step
def fun(x,y,half_life = hl ):
    dydx = np.zeros(len(half_life)+1)
    dydx[0] = -y[0]/half_life[0]
    for i in range(0,len(half_life)-1):
        dydx[i+1] = y[i]/half_life[i] - y[i+1]/half_life[i+1]
    dydx[-1] = y[-2]/half_life[-2]
    return dydx

#setting up the initial conditions of the sets of ODEs
t0 = 0
t1 = 2*10.**21
y0 = [0]*(len(hl)+1)
y0[0] = 1

# setting our timesteps to 1000 points.
t_ev = np.logspace(0, np.log10(0.9*t1), num = 1000)

t_s = time.time()
# ans_rk4 = integrate.solve_ivp(fun,[t0,t1],y0,t_eval=t_ev, method = 'RK45')
# ans_rk2 = integrate.solve_ivp(fun,[t0,t1],y0 ,method = 'RK23')
# both of those methods take too long to run
ans_radau = integrate.solve_ivp(fun, [t0,t1], y0 , t_eval=t_ev, method = 'Radau')

t_en = time.time()


print(t_en - t_s)

# extract the number from the output of the ODE solver
u238 , pb206 = ans_radau.y[0],ans_radau.y[14]
th230, u234= ans_radau.y[4], ans_radau.y[3]
t = ans_radau.t
#
r1 = pb206/u238
r2 = th230/u234

# the last 100 points are of our interest
plt.plot(t[900:],r1[900:],label = 'Pb206/U238')
plt.plot(t[900:],u238[900:],label = 'U238')
plt.plot(t[900:],pb206[900:],'-.', label = 'Pb206')
plt.legend()
plt.xscale('log')
plt.xlabel('Timesteps')
plt.ylabel('Ratio')
plt.savefig("upb.png",dpi=300, bbox_inches = "tight")
plt.show()

# plot the region of interest
plt.plot(t[600:-150],r2[600:-150],label = 'Th230/U234')
#plt.plot(t[600:-150],th230[600:-150])
#plt.plot(t[600:-150], u234[600:-150],'-.')
#plt.yscale('log')
plt.xscale('log')
plt.xlabel('Timesteps')
plt.ylabel('Ratio')
plt.legend()
plt.savefig("thu.png",dpi=300, bbox_inches = "tight")
plt.show()