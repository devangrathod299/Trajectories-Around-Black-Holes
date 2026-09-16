# Trajectories-Around-Black-Holes
Computational study of geodesic motion in curved spacetime, progressing from Schwarzschild trajectories to Kerr (rotating) black holes, followed by photon detection studies from the Schwarzschild ISCO.

## Physics Pipeline:

### 1) Metric

Trajectories are computed in the equatorial plane ($\theta = {\pi}/{2}$). The Schwarzschild metric used:

$$
\begin{aligned}
ds^2 = -(1 - \frac{R}{r}){dt^2} + {(1 - \frac{R}{r})}^{-1}{dr^2} + {r^2}{{d\theta}^{2}} + {r^2}{sin\theta^2}{d\phi^2}
\end{aligned}
$$

where R = $2GM/c^2$ is the Schwarzschild radius. The metric supplies two conserved quantities along any geodesic:

* Angular momentum: L
* Energy: E

For the Kerr extension, the Boyer–Lindquist form of the metric is used, introducing the spin parameter a = J/(Mc) and the frame-dragging term that couples t and $\phi$.

### 2)Geodesic Equation:

The general geodesic equation,

$$
\begin{aligned}
\frac{d^2x^{\alpha}}{d{\lambda}^2} + {\Gamma}^{\alpha}_{\mu\beta} {\left(\frac{dx^{\mu}}{d{\lambda}}\right)}{\left(\frac{dx^{\beta}}{d{\lambda}}\right)} = 0
\end{aligned}
$$

is not integrated directly in Christoffel-symbol form. Instead, E and L are used to reduce the radial motion to a first-order "energy" equation:

**Massless (photon):**

$$
\begin{aligned}
\left(\frac{dr}{d\lambda}\right)^2 = E^2 - \left(1 - \frac{R}{r}\right)\left(\frac{L^2}{r^2}\right)
\end{aligned}
$$

**Massive $(\epsilon = 1):$**

$$
\begin{aligned}
\left(\frac{dr}{d\lambda}\right)^2 = E^2 - \left(1 - \frac{R}{r}\right)\left(\frac{L^2}{r^2}\right) - \epsilon\left(1 - \frac{R}{r}\right)
\end{aligned}
$$

with the azimuthal equation $\frac{d\phi}{d\lambda} = \frac{L}{r^2}$ in both cases.

Since this is a squared first derivative, it is differentiated once more to obtain a clean second-order ODE suitable for RK4:

$$
\begin{aligned}
\frac{d^2r}{d\lambda^2} = \frac{L^2}{r^3} - \left(\frac{3}{2}\frac{RL^2}{r^4}\right) - \epsilon\frac{R}{(2r^2)}
\end{aligned}
$$

For Kerr, the equivalent (E, L)-parametrized radial and azimuthal equations include additional a-dependent terms from frame dragging.

### 3) Initial Conditions

Each trajectory requires:

* $r_0, \phi_0$ — initial radial distance and azimuthal angle (hard-coded)
* ${\frac{dr}{d\lambda}}$ — obtained algebraically from the energy equation:
* $\frac{dr}{d\lambda} = \sqrt{E^2 - \left(1 - \frac{R}{r}\right)\left(\frac{L^2}{r^2}\right) - \epsilon\left(1 - \frac{R}{r}\right)}$

E, L — set by the physical scenario (e.g. circular orbit conditions, or a chosen impact parameter b = L/E for scattering trajectories)

These are fed into an RK4 integrator (later migrated to scipy.solve_ivp) to march the trajectory forward in the affine parameter λ.

### 4) Validation

Trajectories are checked against the weak-field deflection formula:

$$
\begin{aligned}
\Delta\phi_{weak} = \frac{2R}{b} = \frac{4GM}{bc^2}
\end{aligned}
$$

by plotting $\Delta\phi_{sim} / \Delta\phi_{weak}$ vs. b. This ratio → 1 for large b (weak-field limit) and diverges from 1 as b decreases toward the black hole, as expected.

### Results:

* Photon sphere: unstable circular photon orbit at r = 1.5R, confirmed by a constant-r trajectory in the massless case.
* Perihelion shift: bound massive-particle orbits precess per revolution rather than closing, consistent with GR (vs. fixed Newtonian ellipses).
* Kerr prograde/retrograde asymmetry: frame dragging causes prograde photons to linger near the black hole longer than retrograde photons at the same impact parameter; weak-field validation shows the two branches diverging in opposite directions from the non-rotating prediction as b decreases.
* ISCO photon detection: rays emitted from the Schwarzschild ISCO (r = 3R) and traced to a distant observer reproduce the Einstein ring condition $(b > b_{crit} = 3√3 R/2)$, and — once source motion, aberration, and Doppler shift are included — show direction-dependent redshift/blueshift patterns in the detected energy vs. source angle.

### References:

* B. F. Schutz, A First Course in General Relativity
* R. M. Wald, General Relativity, Ch. 12
