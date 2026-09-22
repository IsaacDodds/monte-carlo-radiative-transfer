# Monte Carlo radiative transfer

A C program that follows photons through a plane-parallel scattering slab by Monte Carlo and records the emergent radiation. Coursework from the BSc in Mathematics and Physics, University of Bath.

The program first checks its two samplers for the Rayleigh phase function, rejection sampling and inversion of the cumulative distribution, and writes the samples to `Rejection.txt` and `Cumulative.txt`. It then runs the photon walk, one million photons per run, for a set of optical depths with isotropic and Rayleigh scattering, writes one file per run named `tau_<tau>_<method>.txt`, and prints the fraction of photons absorbed.

## Files

- `radiative-transfer-solver.c`, the source
- `radiative-transfer-solver.pdf`, the report

## Build and run

```
cc -O2 -o radiative radiative-transfer-solver.c -lm
./radiative
```

## Licence

MIT. See `LICENSE`.
