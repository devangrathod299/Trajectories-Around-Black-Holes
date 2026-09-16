"""
Kerr geodesic integrator (Boyer-Lindquist coordinates, equatorial plane).

    d^2r/dlambda^2 = E^2 + a^2*(E^2+epsilon) - L^2/r^2
                      + R*(L-a*E)^2/r^3 - epsilon*(1-R/r)
    dphi/dlambda    = (1/Delta) * (E*R*a/r + L*(1-R/r))

Delta = r^2 - R*r + a^2. Produces prograde/retrograde photon trajectories
and frame-dragging effects (spherical asymmetry, ergosphere behaviour).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Functions 
def f(lam, state, R, L, E, a):
    r, dr_dlam, phi = state
    delta = (r**2) - (R *r) + (a**2)
    return np.array([dr_dlam, ((((L**2) - ((E**2) * (a**2)))/(r**3)) - ((3 * R * ((L - (E * a))**2))/(2 * (r**4)))), ((((E * R * a)/r) + (L * (1 - (R/r))))/delta)])

def make_hit_horizon(R,a):
    def hit_horizion(lam, state):
        r = state[0]
        return r - ((R + np.sqrt((R**2) - (4 * (a**2))))/2)   
    hit_horizion.terminal = True
    hit_horizion.direction = -1   
    return hit_horizion

#Set-up
R, E, a = 1, 1, (-0.5)
L = 3.5

lam_start, lam_end, lam_step = 0, 10000, 100000
lam_range = np.linspace(lam_start, lam_end, lam_step)

r0, phi0 = 5000, 0
dr_dlam_0 = (-np.sqrt((E**2) + (((a**2) - (L**2))/(r0**2)) + ((R * ((L - (a * E))**2))/(r0**3))))

state = np.array([r0, dr_dlam_0, phi0])

#Scipy
f_wrapped = lambda lam, state : f(lam, state, R, L, E, a)
event_fn = make_hit_horizon(R,a)

sol = solve_ivp(f_wrapped, (lam_start, lam_end), state, method = 'DOP853', t_eval = lam_range, dense_output = True, events = (event_fn), atol = 1e-8, rtol = 1e-8)

#Output data accumalation
r_array = sol.y[0]
phi_array = sol.y[2]

#Cartesian coordinates
x = r_array * (np.cos(phi_array))
y = r_array * (np.sin(phi_array))

i_vector_x, i_vector_y = (x[100] - x[0]), (y[100] - y[0])
f_vector_x, f_vector_y = (x[-1] - x[-100]), (y[-1] - y[-100])

mag_i = np.sqrt((i_vector_x**2) + (i_vector_y**2))
mag_f = np.sqrt((f_vector_x**2) + (f_vector_y**2))

iu_vector = np.array([(i_vector_x/mag_i), (i_vector_y/mag_i)])
fu_vector = np.array([(f_vector_x/mag_f), (f_vector_y/mag_f)])

def dot(iu_vector, fu_vector):
    i_f = (iu_vector[0] * fu_vector[0]) + (iu_vector[1] * fu_vector[1])
    deflect = np.degrees(np.arccos(i_f))
    return deflect

delta_phi = np.degrees(((2 * R)/L))

print(f"Simulated deflection - \033[1m{dot(iu_vector, fu_vector)}\033[0m degrees\nweak field - \033[1m{delta_phi}\033[0m degrees.\nratio - \033[1m{dot(iu_vector, fu_vector)/delta_phi}\033[0m")

#Plotting
plt.style.use("dark_background")
fig, axes = plt.subplots(1,2, figsize = (10,4))

axes[0].plot(phi_array, np.log10(r_array), color = "orange")
axes[0].set_xlabel('\u03D5')
axes[0].set_ylabel('log\u2081\u2080r')
axes[0].grid(True)

axes[1].plot(x, y, color = "#00ff73")
axes[1].plot(0,0, color = 'white', marker = 'o')
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')
axes[1].set_xlim(-4,6)
axes[1].set_ylim(-4,4)
axes[1].grid(True)

plt.show()
