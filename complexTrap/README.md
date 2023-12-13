This fold contains the MATLAB code accompanying the paper:

- J. Hoskins, M. Rachh, and B. Wu "On quadrature for singular integral operators with complex symmetric quadratic forms", arXiv preprint arXiv:2312.06812  (2023)

This code contructs corrected trapezoidal rule for the Helmholtz integral operators on "complexified" surfaces in 3D, where the first fundamental form can be copmlex symmetric.

- `test_conv_bump.m` and `test_conv_cyl.m` : convergence tests for corrected trapezoidal rules for Helmholtz layer potentials on a complexified bump and a complexified cylinder. (See Example 1 in the paper.)
- `test_complexTrap.m` : solves a Helmholtz boundary value problem on a Gaussian bump.
- `epstein_zeta_cmpl.m` : evaluation of the Epstein zeta function complex symmetric quadratic forms.

