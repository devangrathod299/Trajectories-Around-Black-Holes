"""
Schwarzschild geodesic integrator.

Integrates the radial ODE for massless (epsilon=0) and massive (epsilon=1)
particles in the equatorial plane, using conserved energy E and angular
momentum L.

    d^2r/dlambda^2 = L^2/r^3 - 1.5*R*L^2/r^4 - epsilon*R/(2*r^2)
    dphi/dlambda    = L/r^2

Initial dr/dlambda is set from:
    dr/dlambda = sqrt(E^2 - (1-R/r)*(L^2/r^2) - epsilon*(1-R/r))

Produces trajectories reproducing the photon sphere (r=1.5R) and
perihelion shift for bound massive orbits.
"""

#Photon trajectory
import numpy as np 
import matplotlib.pyplot as plt

#Functions
def f(a, state, R, L):
    r, v, phi = state
    return np.array([v, ( ((L**2)/(r**3)) - ((3*R*(L**2))/(2*(r**4)))), (L/(r**2))])

def rk4_step(f, a, state, h, R, L):
    k1 = f(a, state, R, L)
    k2 = f((a + (h/2)), (state + ((h/2) * k1)), R, L)
    k3 = f((a + (h/2)), (state + ((h/2) * k2)), R, L)
    k4 = f((a + h), (state + (h * k3)), R, L)

    return (state + ((h * (k1 + (2 * k2) + (2 * k3) + k4)) / 6))

#Set up
#Constants
a_start, a_end,  h = 0, 1000, 0.01
R, L, E = 1, ((3 * np.sqrt(3))/2), 1

#Paramter range
a_values = np.arange(a_start, a_end, h)

#Initial COndition
r0, phi0 = 3, (0)
dr_da_0 = round((-np.sqrt((E**2) - ((1 - (R/r0)) * ((L**2)/(r0**2))))), 4)

state = np.array([r0, dr_da_0, phi0])       

#Looping
r_list = []
phi_list = []

for a in a_values:
    if state[0] < R or state[0] > 1000:
        break
    r_list.append(state[0])
    phi_list.append(state[2])
    state = rk4_step(f, a, state, h, R, L)

#Extra plots
#trajectory
r_array = np.array(r_list)
phi_array = np.array(phi_list)
x = r_array * np.cos(phi_array)
y = r_array * np.sin(phi_array)

#Original trajectory
slope = (y[100]-y[0])/(x[100]-x[0])
p = np.arange(x[-1],x[0],h)
q = (slope * p) + (min(r_list))

angle = np.arange(0, (2 * np.pi), h)
n = np.cos(angle)
m = np.sin(angle)

#After deflection:
x_out1, y_out1 = x[-50], y[-50]
x_out2, y_out2 = x[-1],  y[-1]

slope_out = (y_out2 - y_out1) / (x_out2 - x_out1)
def_x = np.arange(x[-1],x[0],h)
def_y = (slope_out * def_x) + (min(r_array))

#Deflection angle
#delta_phi = round(np.degrees(phi_list[-1] - np.pi), 3)

d_phi = np.degrees(round((np.arctan((slope_out - slope)/(1 + (slope * slope_out)))) , 4))
delta_phi = round(np.degrees((2 * R)/L), 4)

#plotting
plt.style.use("dark_background")

fig, axes = plt.subplots(1,2, figsize = (10,4))

axes[0].plot(phi_array, np.log10(r_array), color = "orange", label = "RK4 simulation plot")
axes[0].set_xlabel('\u03D5')
axes[0].set_ylabel('log\u2081\u2080r')
axes[0].grid(True)
axes[0].legend()

axes[1].plot(x ,y , color = "#2b80ff", label = "Trajecotry")
axes[1].plot(0,0, color = 'white', label = "Center", marker = 'o')
axes[1].plot(p,q, color = "#ffffff", label = 'Original trajectory', linestyle='--', linewidth=0.8)
axes[1].plot(def_x,def_y, color = "#ffffff", label = 'Deflection trajectory', linestyle='--', linewidth=0.8)
axes[1].plot(n,m, color = 'yellow', label = 'Event horizon')
axes[1].set_xlabel('x')
axes[1].set_xlim(-5,5)
axes[1].set_ylabel('y')
axes[1].set_ylim(-5,5)
axes[1].grid(True)
axes[1].legend()

#plt.axis('equal')
plt.show()
print(f"The deviation between the original path and the curved parth is \033[1m{round(d_phi, 4)}\033[0m degrees via RK4 simulation while weak field predicts \033[1m{delta_phi}\033[0m degrees.\n")
