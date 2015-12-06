
#ifndef FCSINTP_H

#include "potential.h"

struct _orbita_models_fcsp3d_trilinear_paramset
{
  // grid points
  float * R_ax,  * Z_ax,  * phi_ax;
  int     R_pts,   Z_pts,   phi_pts;

  // data array
  float * psi_val, * FR_val,
        * FZ_val, * Fphi_val;

  // other params
  float   dphi;
};

double _orbita_models_fcsp3d_trilinear_rho(const void *, const double *);
double _orbita_models_fcsp3d_trilinear_Mr(const void *, const double);
double _orbita_models_fcsp3d_trilinear_Mt(const void *);
double _orbita_models_fcsp3d_trilinear_psi(const void *, const double *);
int _orbita_models_fcsp3d_trilinear_f(const void *, const double *, double *);
int _orbita_models_fcsp3d_kill(struct orbita_potential *);

struct orbita_potential * orbita_potential_fcsp3d_trilinear(const char *, const int, const int, const int, const double, const double, const double);
int _orbita_models_fcsp3d_kill(struct orbita_potential *);

#define FCSINTP_H
#endif
