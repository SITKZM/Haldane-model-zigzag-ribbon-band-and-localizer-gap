# Haldane-model-zigzag-ribbon-band-and-localizer-gap
Calculation of the zigzag ribbon bands of the Haldane model and the localizer gap to distinguish the topological phases.

## Examples
There are three plots as examples, which were created using the following parameters:
Nx = 10
band
    ky = np.linspace(0, 2 * np.pi / np.sqrt(3), 100)
localizer gap, local Chern number
    Ny = 10
    kappa = 1

on-site potential at sublattices a(-), b(+)
m = 0
on-site potential at triangle lattice
m2 = -0.35
nearest hopping
t1 = 1
triangle lattice hopping
t2 = 0.2
hexagonal-triangle hopping
t3 = 0.3
next nearest hopping
tc = 0.5
phi = np.pi / 2


## Reference
A. Cerjan and T. A. Loring. Phys. Rev. B 106, 064109 (2022).
