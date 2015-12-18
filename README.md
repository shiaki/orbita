# orbita

**Numerical library for the study of stellar orbits in galactic potentials.**

(migrated to GitHub in Nov. 2015)

## Author

Yujing Qin, Shanghai Astron. Obs., Chinese Academy of Sciences
contact: qinyj (dot) astro (at) gmail (dot) com

## Introduction

*orbita* is a numerical library (and open framework) to study stellar orbits in galactic potentials. It integrates orbits in various potential models, and performs on-the-fly and post-integration analysis.

### Potential models

*orbita* supports a variety of highly flexible potential models. Potentials can be single- or multi-component, and each component can be analytical, grid-interpolated or in basis-expansion form, with static or time-dependent parameters. Even hierarchically-organized potentials are allowed (but not recommended). Custom potential models can be easily added to the existing framework.

### Integrators

Integrators in *orbita* are also modularized. Currently only 4th order Runge-Kutta integrators are implemented (in both inertia and co-rotating frames), but it is trivial to add other integrators to this framework. Time step can be adaptively adjusted to satisfy the required level of error.

### On-the-fly analyzers

*orbita* has interfaces to perform on-the-fly analysis during integration. Besides doing simple tasks like out-of-bounds detection, analyzers can control orbit integration by writing orbit data and moving the pointer of integration. Many features are implemented as analyzers, including making Poincare sections, searching for periodic orbits, etc.

## Dependencies

[GSL](http://www.gnu.org/software/gsl/) v1.8 or later
Required in: evaluating special functions in some potentials, finding periodic orbits using minimization and root-finding.

## Status

In active development.

Potential models implemented:

* Point mass
* uniform-density sphere
* Plummer sphere
* pseudo-isothermal sphere
* Miyamoto-Nagai disk
* Ferrers n=2 bar
* cylindrical grid interpolator
* spherical grid interpolator

Potential models available in pieces (to be added):

* razor-thin exponential disk
* double exponential disk
* exponential-sech disk
* Gaussian ring profile
* multipole potential

## License

The software is licensed under a 3-clause BSD-style license. See LICENSE.txt for details.
