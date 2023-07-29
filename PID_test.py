"""
Created on Sat Sep 17 23:39:55 2022
https://youtu.be/k46nCvOBllA
http://apmonitor.com/pdc/index.php/Main/FeedbackControl
@author: Nakarin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# define model
def vehicle(v,t,u):
    # inputs
    #  v    = vehicle velocity (m/s)
    #  t    = time (sec)
    #  u    = gas pedal position (-50% to 100%)
    #  load = passenger load + cargo (kg)
    Cd = 0.24    # drag coefficient
    rho = 1.225  # air density (kg/m^3)
    A = 5.0      # cross-sectional area (m^2)
    Fp = 200     # thrust parameter (N/%pedal)
    m = 700      # vehicle mass (kg)
    # calculate derivative of the velocity
    dv_dt = (1.0/(m)) * (Fp*u - 0.5*rho*Cd*A*v**2)
    return dv_dt

tf = 300.0                 # final time for simulation
nsteps = 1000              # number of time steps
delta_t = tf/(nsteps-1)   # how long is each time step?
ts = np.linspace(0,tf,nsteps) # linearly spaced time vector

# simulate step test operation
step = np.zeros(nsteps) # u = valve % open


# velocity initial condition
v0 = 0.0
# set point
sp = 25.0
# for storing the results
vs = np.zeros(nsteps)
sp_store = np.zeros(nsteps)

ubias=0
sum_int=0
us = np.zeros(nsteps)
es = np.zeros(nsteps)
ies = np.zeros(nsteps)

# PID
# IMC Tuning Correlations
k_p=40/25
theta_p=0 # Dead Time
tau_p=24

# Use PI only
tau_c = tau_p
k_c=(1/k_p)*(tau_p/(theta_p+tau_p))
tau_i = tau_p 


#Initial guess for PID
P = k_c
I=k_c/tau_c 


# Adjust PI
P=2
I=0.09

k_c = P
tau_i =k_c/I 

# simulate with ODEINT
for i in range(nsteps-1):
    if i== int(nsteps/5):
        sp = 10.0
    if i== int(2*nsteps/5):
        sp = 40.0
    if i== int(3*nsteps/5):
        sp = 100.0
    if i== int(4*nsteps/5):
        sp = 0.0
    sp_store[i+1] = sp
    error = sp-v0
    es[1+1]=error
    sum_int=sum_int+error*delta_t
    u=ubias+k_c*error+(k_c/tau_i)*sum_int
    
    # Anti wire up
    # clip inputs to -50% to 100%
    if u >= 100.0:
        u = 100.0
    if u <= -50.0:
        u = -50.0
    
    ies[i+1]=sum_int
    es[i+1]=error
    us[i+1]=u
    v = odeint(vehicle,v0,[0,delta_t],args=(u,))
    v0 = v[-1]   # take the last value
    vs[i+1] = v0 # store the velocity for plotting
plt.figure(figsize=(14, 15))
plt.subplot(4,1,1)
plt.plot(ts,vs,'b-',linewidth=3)
plt.plot(ts,sp_store,'k--',linewidth=2)
plt.ylabel('Velocity (m/s)', fontsize=14)
plt.legend(['Velocity','Set Point'],loc=2)
plt.subplot(4,1,2)
plt.plot(ts[0:i+1],us[0:i+1],'r',linewidth=3)
plt.ylabel('Gas Pedal', fontsize=14)    
plt.legend(['Gas Pedal (%)'])

plt.subplot(4,1,3)
plt.plot(ts[0:i+1],es[0:i+1],'b',linewidth=3)
plt.ylabel('Error', fontsize=14)    
plt.legend(['Error'])

plt.subplot(4,1,4)
plt.plot(ts[0:i+1],ies[0:i+1],'g',linewidth=3)
plt.ylabel('Integration', fontsize=14)    
plt.legend(['Integration'])
plt.xlabel('Time (sec)', fontsize=14)
plt.savefig('PI_vehicle.pdf', dpi=150)

