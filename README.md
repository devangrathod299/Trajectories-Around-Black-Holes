# Trajectories-Around-Black-Holes
Computational study of geodesic motion in curved spacetime, progressing from Schwarzschild trajectories to Kerr (rotating) black holes, followed by photon detection studies from the Schwarzschild ISCO.

## Physics Pipeline:

### 1) Metric

Trajectories are computed in the equatorial plane ($\theta = {\pi}/{2}$). The Schwarzschild metric used:

$ds^2 = -(1 - \frac{R}{r}){dt^2} + {(1 - \frac{R}{r})}^{-1}{dr^2} + {r^2}{{d\theta}^{2}} + {r^2}{sin\theta^2}{d\phi^2}$

where R = $2GM/c^2$ is the Schwarzschild radius. The metric supplies two conserved quantities along any geodesic:
