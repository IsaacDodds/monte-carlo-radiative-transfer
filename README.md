# Radiative Transfer Solver

This repository contains a C-based implementation of a 1D radiative transfer solver. It simulates photon transport through a scattering medium using both deterministic PDE discretisation and Monte Carlo methods.

## Features

- Solves the time-independent radiative transfer equation
- Isotropic and Rayleigh scattering models
- Finite-difference discretisation
- Monte Carlo rejection and inverse transform sampling
- Emergent intensity and angular flux output
- Adjustable opacity, density, and slab width

## Files

- `id408.c` — Main C source code
- `id408.pdf` — Report summarising physics and results
- `*.txt` — Output intensity or flux files (optional)

## Author

Isaac Dodds  
Physics + Mathematics | University of Bath  
Commonwealth Games athlete and scientific computing student  
[GitHub Profile](https://github.com/IsaacDodds)

## License

MIT License
