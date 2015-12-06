
#ifndef SPHINTP_H

struct _orbita_models_sphintp_trilinear_paramset
{
  // grid points
  double * r_ax,  * theta_ax,  * phi_ax;
  int      r_pts,   theta_pts,   phi_pts;

  // data array
  double * psi_val, * Fx_val, * Fy_val,  * Fz_val;

  // other params
  double   dphi;
};

double  _orbita_models_sphintp_trilinear_rho(const void *, const double *);
double  _orbita_models_sphintp_trilinear_Mr(const void *, const double);
double  _orbita_models_sphintp_trilinear_Mt(const void *);
double  _orbita_models_sphintp_trilinear_psi(const void *, const double *);
int     _orbita_models_sphintp_trilinear_f(const void *, const double *, double *);
int     _orbita_models_sphintp_kill(struct orbita_potential *);

struct orbita_potential *
orbita_potential_sphintp_trilinear(const char *, const int, const int, const int, const double);
int     _orbita_models_sphintp_kill(struct orbita_potential *);

#define SPHINTP_H
#endif
