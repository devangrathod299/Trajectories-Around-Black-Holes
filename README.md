# Trajectories-Around-Black-Holes
Computational study of geodesic motion in curved spacetime, progressing from Schwarzschild trajectories to Kerr (rotating) black holes, followed by photon detection studies from the Schwarzschild ISCO.

## Physics Pipeline:

### 1) Metric

Trajectories are computed in the equatorial plane ($\theta = {\pi}/{2}$). The Schwarzschild metric used:

$ds^2 = -(1 - \frac{R}{r}){dt^2} + {(1 - \frac{R}{r})}^{-1}{dr^2} + {r^2}{{d\theta}^{2}} + {r^2}{sin\theta^2}{d\phi^2}$

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
\left(\frac{dr}{d\lambda}\right)^2 = E^2 - \left(1 - \frac{R}{r})(\frac{L^2}{r^2}\right)
\end{aligned}
$$

**Massive $(\epsilon = 1):$**

$$
\begin{aligned}
\left(\frac{dr}{d\lambda}\right)^2 = E^2 - \left(1 - \frac{R}{r}\right)\left(\frac{L^2}{r^2}\right) - \epsilon\left(1 - \frac{R}{r}\right)$

with the azimuthal equation $\frac{d\phi}{d\lambda} = \frac{L}{r^2} in both cases.

Since this is a squared first derivative, it is differentiated once more to obtain a clean second-order ODE suitable for RK4:

$$
\begin{aligned}
\frac{d^2r}{d\lambda^2} = \frac{L^2}{r^3} - (\frac{3}{2}) R \frac{L^2}{r^4} - \epsilon\frac{R}{(2r^2)}
\end{aligned}
$$

For Kerr, the equivalent (E, L)-parametrized radial and azimuthal equations include additional a-dependent terms from frame dragging.
