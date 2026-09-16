"""
Compares prograde vs retrograde Kerr photon trajectories at fixed impact
parameter b, and validates both against the (rotation-blind) weak-field
prediction, showing the two branches diverging in opposite directions as
b decreases (frame-dragging asymmetry).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Functions 
def f_a(lam, state, R, L, E, a):
    r, dr_dlam, phi = state
    delta = (r**2) - (R *r) + (a**2)
    return np.array([dr_dlam, ((((L**2) - ((E**2) * (a**2)))/(r**3)) - ((3 * R * ((L - (E * a))**2))/(2 * (r**4)))), ((((E * R * a)/r) + (L * (1 - (R/r))))/delta)])

def make_hit_horizon_a(R,a):
    def hit_horizion_a(lam, state):
        r = state[0]
        return r - ((R + np.sqrt((R**2) - (4 * (a**2))))/2)   
    hit_horizion_a.terminal = True
    hit_horizion_a.direction = -1   
    return hit_horizion_a

def deflection_a(L, R, E):

    lam_start, lam_end, lam_step = 0, 10000, 100000
    lam_range = np.linspace(lam_start, lam_end, lam_step)

    r0, phi0 = 5000, 0
    dr_dlam_0 = (-np.sqrt((E**2) + (((a**2) - (L**2))/(r0**2)) + ((R * ((L - (a * E))**2))/(r0**3))))

    state = np.array([r0, dr_dlam_0, phi0])

    #Scipy
    f_wrapped = lambda lam, state : f_a(lam, state, R, L, E, a)
    event_fn = make_hit_horizon_a(R,a)

    sol_a = solve_ivp(f_wrapped, (lam_start, lam_end), state, method = 'DOP853', t_eval = lam_range, dense_output = True, events = (event_fn), atol = 1e-8, rtol = 1e-8)

    #Output data accumalation
    r_array = sol_a.y[0]
    phi_array = sol_a.y[2]

    #Cartesian coordinates
    x = r_array * (np.cos(phi_array))
    y = r_array * (np.sin(phi_array))

    i_vector_x, i_vector_y = (x[100] - x[0]), (y[100] - y[0])
    f_vector_x, f_vector_y = (x[-1] - x[-100]), (y[-1] - y[-100])

    mag_i = np.sqrt((i_vector_x**2) + (i_vector_y**2))
    mag_f = np.sqrt((f_vector_x**2) + (f_vector_y**2))

    iu_vector = np.array([(i_vector_x/mag_i), (i_vector_y/mag_i)])
    fu_vector = np.array([(f_vector_x/mag_f), (f_vector_y/mag_f)])

    i_f = (iu_vector[0] * fu_vector[0]) + (iu_vector[1] * fu_vector[1])
    def_lo = np.degrees(np.arccos(i_f))

    def_wf = np.degrees(((2 * R)/L))
    ratio = def_lo/def_wf

    return ratio

def f_b(lam, state, R, L, E, b):
    r, dr_dlam, phi = state
    delta = (r**2) - (R *r) + (b**2)
    return np.array([dr_dlam, ((((L**2) - ((E**2) * (b**2)))/(r**3)) - ((3 * R * ((L - (E * b))**2))/(2 * (r**4)))), ((((E * R * b)/r) + (L * (1 - (R/r))))/delta)])

def make_hit_horizon_b(R,b):
    def hit_horizion_b(lam, state):
        r = state[0]
        return r - ((R + np.sqrt((R**2) - (4 * (b**2))))/2)   
    hit_horizion_b.terminal = True
    hit_horizion_b.direction = -1   
    return hit_horizion_b

def deflection_b(L, R, E):

    lam_start, lam_end, lam_step = 0, 10000, 100000
    lam_range = np.linspace(lam_start, lam_end, lam_step)

    r0, phi0 = 5000, 0
    dr_dlam_0 = (-np.sqrt((E**2) + (((b**2) - (L**2))/(r0**2)) + ((R * ((L - (b * E))**2))/(r0**3))))

    state = np.array([r0, dr_dlam_0, phi0])

    #Scipy
    f_wrapped = lambda lam, state : f_b(lam, state, R, L, E, b)
    event_fn = make_hit_horizon_b(R,b)

    sol_b = solve_ivp(f_wrapped, (lam_start, lam_end), state, method = 'DOP853', t_eval = lam_range, dense_output = True, events = (event_fn), atol = 1e-8, rtol = 1e-8)

    #Output data accumalation
    r_array = sol_b.y[0]
    phi_array = sol_b.y[2]

    #Cartesian coordinates
    x = r_array * (np.cos(phi_array))
    y = r_array * (np.sin(phi_array))

    i_vector_x, i_vector_y = (x[100] - x[0]), (y[100] - y[0])
    f_vector_x, f_vector_y = (x[-1] - x[-100]), (y[-1] - y[-100])

    mag_i = np.sqrt((i_vector_x**2) + (i_vector_y**2))
    mag_f = np.sqrt((f_vector_x**2) + (f_vector_y**2))

    iu_vector = np.array([(i_vector_x/mag_i), (i_vector_y/mag_i)])
    fu_vector = np.array([(f_vector_x/mag_f), (f_vector_y/mag_f)])

    i_f = (iu_vector[0] * fu_vector[0]) + (iu_vector[1] * fu_vector[1])
    def_lo = np.degrees(np.arccos(i_f))

    def_wf = np.degrees(((2 * R)/L))
    ratio = def_lo/def_wf

    return ratio

def f_c(lam, state, R, L, E, c):
    r, dr_dlam, phi = state
    delta = (r**2) - (R *r) + (c**2)
    return np.array([dr_dlam, ((((L**2) - ((E**2) * (c**2)))/(r**3)) - ((3 * R * ((L - (E * c))**2))/(2 * (r**4)))), ((((E * R * c)/r) + (L * (1 - (R/r))))/delta)])

def make_hit_horizon_c(R,c):
    def hit_horizion_c(lam, state):
        r = state[0]
        return r - ((R + np.sqrt((R**2) - (4 * (c**2))))/2)   
    hit_horizion_c.terminal = True
    hit_horizion_c.direction = -1   
    return hit_horizion_c

def deflection_c(L, R, E):

    lam_start, lam_end, lam_step = 0, 10000, 100000
    lam_range = np.linspace(lam_start, lam_end, lam_step)

    r0, phi0 = 5000, 0
    dr_dlam_0 = (-np.sqrt((E**2) + (((c**2) - (L**2))/(r0**2)) + ((R * ((L - (c * E))**2))/(r0**3))))

    state = np.array([r0, dr_dlam_0, phi0])

    #Scipy
    f_wrapped = lambda lam, state : f_c(lam, state, R, L, E, c)
    event_fn = make_hit_horizon_c(R,c)

    sol_c = solve_ivp(f_wrapped, (lam_start, lam_end), state, method = 'DOP853', t_eval = lam_range, dense_output = True, events = (event_fn), atol = 1e-8, rtol = 1e-8)

    #Output data accumalation
    r_array = sol_c.y[0]
    phi_array = sol_c.y[2]

    #Cartesian coordinates
    x = r_array * (np.cos(phi_array))
    y = r_array * (np.sin(phi_array))

    i_vector_x, i_vector_y = (x[100] - x[0]), (y[100] - y[0])
    f_vector_x, f_vector_y = (x[-1] - x[-100]), (y[-1] - y[-100])

    mag_i = np.sqrt((i_vector_x**2) + (i_vector_y**2))
    mag_f = np.sqrt((f_vector_x**2) + (f_vector_y**2))

    iu_vector = np.array([(i_vector_x/mag_i), (i_vector_y/mag_i)])
    fu_vector = np.array([(f_vector_x/mag_f), (f_vector_y/mag_f)])

    i_f = (iu_vector[0] * fu_vector[0]) + (iu_vector[1] * fu_vector[1])
    def_lo = np.degrees(np.arccos(i_f))

    def_wf = np.degrees(((2 * R)/L))
    ratio = def_lo/def_wf

    return ratio

R, E, a, b, c = 1, 1, (-0.5), 0, 0.5

L_start, L_end, step = 4, 1000, 1
L_range = np.arange(L_start, L_end, step)

L_list_a = []
for L in L_range:
    L_list_a.append(deflection_a(L, R, E))
    L = L + 1

L_list_b = []
for L in L_range:
    L_list_b.append(deflection_b(L, R, E))
    L = L + 1

L_list_c = []
for L in L_range:
    L_list_c.append(deflection_c(L, R, E))
    L = L + 1

#Plotting
plt.style.use("dark_background")
plt.plot(np.log10(L_range), L_list_a, color = 'orange', label = 'Prograde')
plt.plot(np.log10(L_range), L_list_b, color = 'pink', label = 'Schwarzschild like')
plt.plot(np.log10(L_range), L_list_c, color = 'Yellow', label = 'Retrograde')
plt.xlabel('log\u2081\u2080(b/R)')
plt.ylabel("\u0394\u03a6 - (Simulated/Weak-field app.)")
plt.grid(True)
plt.legend(loc = 'upper right')
plt.show()
