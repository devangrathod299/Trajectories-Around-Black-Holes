"""
Extends forward ray-tracing to a moving source (circular ISCO orbit).
Applies special-relativistic aberration and Doppler shift on top of
gravitational redshift to compute observed photon energy vs. source
angular position.

    v_obs = v_emit * ((1-beta*cos(theta_o))/sqrt(1-beta^2))
            * sqrt((1-R/r0)/(1-R/r_obs))

Tangential velocity at ISCO: beta = 1/2 (Schwarzschild ISCO result).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#functions 
def f(lam, state, R, L):
    r, dr_dlam, phi = state
    return np.array([dr_dlam,(((L**2)/(r**3)) - ((3 * R * (L**2))/(2 * (r**4)))), (L/(r**2))])

def plot_neg_18(b, R, E, r0, phi0):
    #se-up
    L = b * E

    dr_dlam_0 = (-np.sqrt((E**2) - (((L**2)/(r0**2)) * (1 - (R/r0))))) 
    
    state = ([r0, dr_dlam_0, phi0])

    lam_start, lam_end, h = 0, 110, 0.001
    lam_range = np.arange(lam_start, lam_end, h)

    f_wrapped = lambda lam, state : f(lam, state, R, L)

    sol = solve_ivp(f_wrapped, (lam_start,lam_end), state, method = 'LSODA', t_eval = lam_range, events = (hit_horizon), rtol = 1e-8, atol = 1e-8)

    r_array = sol.y[0]
    phi_array = sol.y[2]

    x = r_array * (np.cos(phi_array))
    y = r_array * (np.sin(phi_array))

    return np.array([x, y])

def plot_pos_18(b, R, E, r0, phi0):
    #se-up
    L = b * E

    dr_dlam_0 = (np.sqrt((E**2) - (((L**2)/(r0**2)) * (1 - (R/r0)))))
    
    state = ([r0, dr_dlam_0, phi0])

    lam_start, lam_end, h = 0, 110, 0.01
    lam_range = np.arange(lam_start, lam_end, h)

    f_wrapped = lambda lam, state : f(lam, state, R, L)

    sol = solve_ivp(f_wrapped, (lam_start,lam_end), state, method = 'LSODA', t_eval = lam_range, events = (hit_horizon), rtol = 1e-8, atol = 1e-8)

    r_array = sol.y[0]
    phi_array = sol.y[2]

    x = r_array * (np.cos(phi_array))
    y = r_array * (np.sin(phi_array))

    return np.array([x, y])

def hit_horizon(lam, state):
    r = state[0]
    return r - R   # zero when r = R

hit_horizon.terminal = True
hit_horizon.direction = -1   

#For light
R,E, = 1,1

#For source
r0 = (3 * R)
r_metric = 1/(np.sqrt(1 - (R/r0)))

phi0_start, phi0_end, phi0_step = 0, (2 * np.pi), (np.pi/18)
phi0_range = np.arange(phi0_start, phi0_end, phi0_step)

E_isco, L_isco = ((2 * np.sqrt(2))/3), (np.sqrt(3))
beta = -((L_isco * np.sqrt(1 - (R/r0)))/(E_isco * r0))
gamma = 1/(np.sqrt(1 - (beta**2)))

#Lists
alpha_start, alpha_end, alpha_step = (-np.pi) , (0), (np.pi/36)
alpha_source = np.arange(alpha_start, alpha_end, alpha_step) 

cos_alpha = (np.cos(alpha_source) - beta)/(1 - (beta * np.cos(alpha_source)))
sin_alpha = (np.sin(alpha_source)/(gamma * (1 - (beta * np.cos(alpha_source)))))
alpha_real =  np.arctan2(sin_alpha, cos_alpha) + (np.pi/2)

#Looping
b_obs, phi0_obs = [], []

for phi0 in phi0_range:

    for alpha in alpha_real:
        b = (r0 * (np.sin(alpha)))*(r_metric)

        full_sol_n = plot_neg_18(b, R, E, r0, phi0)
        x_n, y_n = full_sol_n[0], full_sol_n[1]

        mask = (x_n >= -115) & (x_n <= -85) & (y_n >= -15) & (y_n <= 15)

        if np.any(mask):
            b_obs.append(b)
            phi0_obs.append(phi0)

        full_sol_p = plot_pos_18(b, R, E, r0, phi0)
        x_p, y_p = full_sol_p[0], full_sol_p[1]

        mask = (x_p >= -110) & (x_p <= -90) & (y_p >= -10) & (y_p <= 10)

        if np.any(mask):
            b_obs.append(b)
            phi0_obs.append(phi0)
        
#Plotting
print(np.degrees(phi0_obs))
plt.style.use("dark_background")
plt.plot(np.degrees(phi0_obs), b_obs, 'o')
plt.xlabel('phi\u2080 of source')
plt.ylabel('impact parameter b of received ray')
plt.grid(True)
plt.show()
