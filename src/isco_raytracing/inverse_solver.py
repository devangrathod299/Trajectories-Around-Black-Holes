"""
Inverse ray-tracing solver: given an observer position (r_obs, phi_obs),
find the impact parameter(s) b of photon(s) reaching that position.

Two-stage pipeline:
  1. Interpolate phi_hit at r_obs for a family of trajectories parametrized by b.
  2. Invert the (b, phi_hit) map to find b at a target phi_obs.

Note: Delta_phi(b) diverges logarithmically as b -> b_crit, producing a
multi-branch structure (photon winding). Each additional 2*pi of winding
requires |b - b_crit| to shrink by ~e^(2*pi) ~ 535, so b sampling near
b_crit must be log-spaced across many decades, not linear.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

#functions 
def f(lam, state, R, L):
    r, dr_dlam, phi = state
    return np.array([dr_dlam,(((L**2)/(r**3)) - ((3 * R * (L**2))/(2 * (r**4)))), (L/(r**2))])

def negative(b, R, E, r0, phi0):
    #set-up
    L = b * E

    dr_dlam_0 = (-np.sqrt((E**2) - (((L**2)/(r0**2)) * (1 - (R/r0))))) 
    
    state = ([r0, dr_dlam_0, phi0])

    lam_start, lam_end, h = 0, 2000, 0.001
    lam_range = np.arange(lam_start, lam_end, h)

    f_wrapped = lambda lam, state : f(lam, state, R, L)

    sol = solve_ivp(f_wrapped, (lam_start,lam_end), state, method = 'LSODA', t_eval = lam_range, events = hit_horizon, rtol = 1e-13, atol = 1e-13)

    r_array = sol.y[0]
    phi_array = sol.y[2]

    return np.array([r_array, phi_array])

def positive(b, R, E, r0, phi0):
    #set-up
    L = b * E

    dr_dlam_0 = (np.sqrt((E**2) - (((L**2)/(r0**2)) * (1 - (R/r0)))))
    
    state = ([r0, dr_dlam_0, phi0])

    lam_start, lam_end, h = 0, 2000, 0.01
    lam_range = np.arange(lam_start, lam_end, h)

    f_wrapped = lambda lam, state : f(lam, state, R, L)

    sol = solve_ivp(f_wrapped, (lam_start,lam_end), state, method = 'LSODA', t_eval = lam_range, events = hit_horizon, rtol = 1e-13, atol = 1e-13)

    r_array = sol.y[0]
    phi_array = sol.y[2]

    return np.array([r_array, phi_array])

def hit_horizon(lam, state):
    r = state[0]
    return r - R   # zero when r = R

hit_horizon.terminal = True
hit_horizon.direction = -1   

#Set-up: Variables
R, E = 1, 1
r0, phi0 = (3 * R), ((0 * np.pi) / 18)

r_obs, phi_obs = 100, ((16.7777 * np.pi)/18)
r_metric = 1/(np.sqrt(1 - (R/r0)))

E_isco, L_isco = ((2 * np.sqrt(2))/3), (np.sqrt(3))
beta = -((L_isco * np.sqrt(1 - (R/r0)))/(E_isco * r0))
gamma = 1/(np.sqrt(1 - (beta**2)))

#Aberration
alpha_start, alpha_end, alpha_step = (-np.pi) , ((1 * np.pi)/18), (np.pi/36)
alpha_source = np.arange(alpha_start, alpha_end, alpha_step) 

alpha_source_crit = np.arcsin(1/np.sqrt(3))

cos_alpha = (np.cos(alpha_source) - beta)/(1 - (beta * np.cos(alpha_source)))
sin_alpha = (np.sin(alpha_source)/(gamma * (1 - (beta * np.cos(alpha_source)))))
alpha_real =  np.arctan2(sin_alpha, cos_alpha) + (np.pi/2)

b_list = []
for alpha in alpha_real:
    b = (r0 * (np.sin(alpha)))*(r_metric)
    b_list.append(b)

#Filtering:
phi_filter, valid_b = [], []

for b in b_list: 
    neg = negative(b, R, E, r0, phi0)
    r_neg, phi_neg = neg[0], neg[1]

    if (r_obs < r_neg.min() or r_obs > r_neg.max()):
        continue 
    
    neg_interpolate = interp1d(r_neg, phi_neg)
    phi_hit_neg = neg_interpolate(r_obs)
    
    phi_filter.append(phi_hit_neg)
    valid_b.append(b)
    

for b in b_list: 
    pos = positive(b, R, E, r0, phi0)
    r_pos, phi_pos = pos[0], pos[1]

    if (r_obs < r_pos.min() or r_obs > r_pos.max()):
        continue  
    
    pos_interpolate = interp1d(r_pos, phi_pos)
    phi_hit_pos = pos_interpolate(r_obs)
    
    phi_filter.append(phi_hit_pos)
    valid_b.append(b)


g = interp1d(phi_filter, valid_b)

print(min(b_list), max(b_list), len(b_list))

print(f"phi filter - {min(phi_filter), max(phi_filter), len(phi_filter)}")

b_sol_pos = g(phi_obs)
b_sol_neg = g(-((2 * np.pi) - phi_obs))

print(f"Solution - {b_sol_pos}")

print(f"alpha - {np.degrees(min(alpha_real)), np.degrees(max(alpha_real))}")

plt.style.use('dark_background')
plt.scatter(np.degrees(alpha_real), b_list)
plt.scatter(np.degrees(phi_filter), valid_b)
plt.grid(True)
plt.show()
