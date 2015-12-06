
#ifndef ORBITA_MODELS_H

#include "potential.h"

/* Many mass models */

/*
  00. Point Mass
    Parameters:
      Mass
*/

double _orbita_models_pointmass_rho(const void *, const double *);
double _orbita_models_pointmass_psi(const void *, const double *);
int    _orbita_models_pointmass_f  (const void *, const double *, double *);
double _orbita_models_pointmass_Mr (const void *, double);
double _orbita_models_pointmass_Mt (const void *);
struct orbita_potential * orbita_potential_pointmass(double);
struct orbita_potential * orbita_potential_pointmass_var(double,
    int (*)(void const * const, const double, void *));
int    _orbita_models_pointmass_kill(struct orbita_potential *);

/* Homogeneous Sphere */

/*
  01. Homogeneous sphere
    Parameters:
      central density rho_0
      radius
*/

double _orbita_models_homogeneous_sphere_rho(const void *, const double *);
double _orbita_models_homogeneous_sphere_psi(const void *, const double *);
int    _orbita_models_homogeneous_sphere_f  (const void *, const double *, double *);
double _orbita_models_homogeneous_sphere_Mr (const void *, double);
double _orbita_models_homogeneous_sphere_Mt (const void *);
struct orbita_potential * orbita_potential_homogeneous_sphere(double, double);
struct orbita_potential * orbita_potential_homogeneous_sphere_var(double, double,
    int (*)(void const * const, const double, void *));
int    _orbita_models_homogeneous_sphere_kill(struct orbita_potential *);

/* Plummer */

/*
  02. Homogeneous sphere
    Parameters:
      total mass
      softening radius
*/

double _orbita_models_plummer_rho(const void *, const double *);
double _orbita_models_plummer_psi(const void *, const double *);
int    _orbita_models_plummer_f  (const void *, const double *, double *);
double _orbita_models_plummer_Mr (const void *, double);
double _orbita_models_plummer_Mt (const void *);
struct orbita_potential * orbita_potential_plummer(double, double);
struct orbita_potential * orbita_potential_plummer_var(double, double,
    int (*)(void const * const, const double, void *));
int    _orbita_models_plummer_kill(struct orbita_potential *);

/* pseudo-isothermal sphere */

double _orbita_models_pseudo_isothermal_rho(const void *, const double *);
double _orbita_models_pseudo_isothermal_psi(const void *, const double *);
int    _orbita_models_pseudo_isothermal_f  (const void *, const double *, double *);
double _orbita_models_pseudo_isothermal_Mr (const void *, double);
double _orbita_models_pseudo_isothermal_Mt (const void *);
struct orbita_potential * orbita_potential_pseudo_isothermal(double, double);
int    _orbita_models_pseudo_isothermal_kill(struct orbita_potential *);

#define ORBITA_MODELS_H
#endif
