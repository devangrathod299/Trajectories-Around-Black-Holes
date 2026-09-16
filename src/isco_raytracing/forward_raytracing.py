"""
Forward ray-tracing from a stationary source at the Schwarzschild ISCO
(r=3R) to a distant observer. Emits N rays at equal angular intervals,
detects which fall within a circular region around the observer, and
records their impact parameters b.

Critical impact parameter for image formation:
    b_crit = 3*sqrt(3)*R/2
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
  
R,E = 1,1

r0 = (3 * R)
phi0_start, phi0_end, phi0_step = 0, (2 * np.pi), (np.pi/18)
phi0_range = np.arange(phi0_start, phi0_end, phi0_step)

theta_start, theta_end, theta_step = (-np.pi/2) , ((5 * np.pi)/9), (np.pi/36)
theta_range = np.arange(theta_start, theta_end, theta_step)

b_obs, phi0_obs = [], []

for phi0 in phi0_range:

    for theta in theta_range:
        b = (r0 * (np.sin(theta)))*(np.sqrt(3/2))

        full_sol_n = plot_neg_18(b, R, E, r0, phi0)
        x_n, y_n = full_sol_n[0], full_sol_n[1]

        mask = (x_n >= -102.5) & (x_n <= -92.5) & (y_n >= -2.5) & (y_n <= 2.5)

        if np.any(mask):
            b_obs.append(b)
            phi0_obs.append(phi0)

        full_sol_p = plot_pos_18(b, R, E, r0, phi0)
        x_p, y_p = full_sol_p[0], full_sol_p[1]

        mask = (x_p >= -105) & (x_p <= -95) & (y_p >= -5) & (y_p <= 5)

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

print(b_obs)
