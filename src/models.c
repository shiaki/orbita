
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "potential.h"

#define EPS 1E-3
#define Sq(X) ((X) * (X))
#define Cb(X) ((X) * (X) * (X))

#define SphRad(R)   (sqrt((R)[0]*(R)[0]+(R)[1]*(R)[1]+(R)[2]*(R)[2]))
#define SphRadSq(R) ((R)[0]*(R)[0]+(R)[1]*(R)[1]+(R)[2]*(R)[2])
#define CylRad(R)   (sqrt((R)[0]*(R)[0]+(R)[1]*(R)[1]))

//const double Pi = 3.14159265358979;
//#define Pi (3.14159265358979)

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
/* 00. Point Mass                                                            */

double
_orbita_models_pointmass_rho(const void * M,
                             const double * R)
{
  return 0.;
}

double
_orbita_models_pointmass_psi(const void * M,
                             const double * R)
{
  double Rt = sqrt(Sq(R[0]) + Sq(R[1]) + Sq(R[2]));
  return -*((double *) M) / Rt;
}

int
_orbita_models_pointmass_f(const void * M,
                           const double * R,
                           double * F)
{
  double Rsq = Sq(R[0]) + Sq(R[1]) + Sq(R[2]);
  double Ft_Rt = *((double *) M) / Rsq;

  if(! isfinite(Ft_Rt)) Ft_Rt = 0;

  F[0] = -Ft_Rt * R[0],
  F[1] = -Ft_Rt * R[1],
  F[2] = -Ft_Rt * R[2];

  return 0;
}

double
_orbita_models_pointmass_Mr(const void * M,
                            double R)
{
  return *((double *) M);
}

double
_orbita_models_pointmass_Mt(const void * M)
{
  return *((double *) M);
}

int
_orbita_models_pointmass_kill(struct orbita_potential * psi)
{
  if(psi -> is_variable)
    free(psi -> init_param);
  free(psi -> param);
  return 0;
}

struct orbita_potential *
orbita_potential_pointmass(double M)
{
  struct orbita_potential * psi =
      (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> param                = (void *) malloc(sizeof(double));
  *((double *)(psi -> param)) = M;

  psi -> symmetry    = ORBITA_SYMM_SPHERICAL;
  psi -> coordinates = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_models_pointmass_rho;
  psi -> psi = & _orbita_models_pointmass_psi;
  psi -> f   = & _orbita_models_pointmass_f;
  psi -> Mr  = & _orbita_models_pointmass_Mr;
  psi -> Mt  = & _orbita_models_pointmass_Mt;

  psi -> is_variable = 0;
  psi -> init_param  = NULL;
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_models_pointmass_kill;

  return psi;
}

struct orbita_potential *
orbita_potential_pointmass_var(double  M,
                               int  (* Param_func)(void const * const,    // p(0), init param
                                                   const double,          // t   , time
                                                   void *))               // p(t), param at time t
{
  struct orbita_potential * psi =
      (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> param                     = (void *) malloc(sizeof(double));
  psi -> init_param                = (void *) malloc(sizeof(double));
  *((double *)(psi -> init_param)) = M;

  psi -> symmetry    = ORBITA_SYMM_SPHERICAL;
  psi -> coordinates = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_models_pointmass_rho;
  psi -> psi = & _orbita_models_pointmass_psi;
  psi -> f   = & _orbita_models_pointmass_f;
  psi -> Mr  = & _orbita_models_pointmass_Mr;
  psi -> Mt  = & _orbita_models_pointmass_Mt;

  psi -> is_variable = 1;
  psi -> param_func  = Param_func;

  return psi;
}

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
/* 01. homogeneous sphere                                                    */

double
_orbita_models_homogeneous_sphere_rho(const void * par,
                                      const double * R)
{
  double Rt = sqrt(Sq(R[0]) + Sq(R[1]) + Sq(R[2]));
  if(Rt < (((double *)par)[1])) return ((double *)par)[0];
  else return 0.;
}

double
_orbita_models_homogeneous_sphere_psi(const void * par,
                                      const double * R)
{
  double Rt = SphRad(R), a = ((double *)par)[1], rho_0 = ((double *)par)[0];
  if(Rt < a) return -2. * Pi * rho_0 * (Sq(a) - (1. / 3.) * SphRadSq(R));
  else return -4. * Pi * rho_0 * Cb(a) / (3. * Rt);
}

double
_orbita_models_homogeneous_sphere_Mr(const void * par,
                                     double R)
{
  double a = ((double *)par)[1], rho_0 = ((double *)par)[0];
  if(R < a) return (4. / 3.) * Pi * rho_0 * Cb(R);
  else return (4. / 3.) * Pi * rho_0 * Cb(a);
}

double
_orbita_models_homogeneous_sphere_Mt(const void * par)
{
  double a = ((double *)par)[1], rho_0 = ((double *)par)[0];
  return (4. / 3.) * Pi * rho_0 * Cb(a);
}

int
_orbita_models_homogeneous_sphere_f(const void * par,
                                    const double * R,
                                    double * F)
{
  double Rsq = Sq(R[0]) + Sq(R[1]) + Sq(R[2]),
         Rt  = sqrt(Rsq);

  double Fr_Rt = -_orbita_models_homogeneous_sphere_Mr(par, Rt) / Rsq / Rt;
  if(! isfinite(Fr_Rt)) Fr_Rt = 0.;

  F[0] = Fr_Rt * R[0],
  F[1] = Fr_Rt * R[1],
  F[2] = Fr_Rt * R[2];

  return 0;
}

int
_orbita_models_homogeneous_sphere_kill(struct orbita_potential * psi)
{
  if(psi -> is_variable)
    free(psi -> init_param);
  free(psi -> param);
  return 0;
}

struct orbita_potential *
orbita_potential_homogeneous_sphere(double rho,
                                    double a)
{
  struct orbita_potential * psi =
      (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> param = (void *) malloc(2 * sizeof(double));
  *((double *)(psi -> param)) = rho;
  *((double *)(psi -> param) + 1) = a;

  psi -> symmetry    = ORBITA_SYMM_SPHERICAL;
  psi -> coordinates = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_models_homogeneous_sphere_rho;
  psi -> psi = & _orbita_models_homogeneous_sphere_psi;
  psi -> f   = & _orbita_models_homogeneous_sphere_f;
  psi -> Mr  = & _orbita_models_homogeneous_sphere_Mr;
  psi -> Mt  = & _orbita_models_homogeneous_sphere_Mt;

  psi -> is_variable = 0;
  psi -> init_param  = NULL;
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_models_homogeneous_sphere_kill;

  return psi;
}

struct orbita_potential *
orbita_potential_homogeneous_sphere_var(double rho,
                                        double a)
{
  orbita_err("orbita_potential_homogeneous_sphere_var: Not implemented.");
  return NULL;
}

/* Plummer */

double _orbita_models_plummer_rho(const void * par, const double * pos)
{
  double M = *((double *)par), b = *((double *)par + 1);
  double R_sq = SphRadSq(pos);
  return (3 * M) * pow(1. + R_sq / Sq(b), -2.5) / (4 * Pi * Cb(b));
}

double _orbita_models_plummer_psi(const void * par, const double * pos)
{
  double M = *((double *)par), b = *((double *)par + 1);
  double R_sq = SphRadSq(pos);
  return -M / sqrt(R_sq + Sq(b));
}

int _orbita_models_plummer_f(const void * par, const double * pos, double * F)
{
  double M = *((double *)par), b = *((double *)par + 1);
  double R_sq = SphRadSq(pos); 
  double Fr_R = -M / pow(Sq(b) + R_sq, 1.5);

  if(! isfinite(Fr_R)) Fr_R = 0;

  F[0] = Fr_R * pos[0],
  F[1] = Fr_R * pos[1],
  F[2] = Fr_R * pos[2];

  return 0;
}

double _orbita_models_plummer_Mr(const void * par, const double R)
{
  double M = *((double *)par), b = *((double *)par + 1);
  return M * Sq(R) * R / pow(Sq(b) + Sq(R), 1.5);
}

double _orbita_models_plummer_Mt(const void * par)
{
  double M = *((double *)par);
  return M;
}

int
_orbita_models_plummer_kill(struct orbita_potential * psi)
{
  if(psi -> is_variable)
    free(psi -> init_param);
  free(psi -> param);
  return 0;
}

struct orbita_potential * orbita_potential_plummer(double M, double b)
{
  struct orbita_potential * psi =
    (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> param = (void *) malloc(2 * sizeof(double));
  *((double *)(psi -> param))     = M;
  *((double *)(psi -> param) + 1) = b;

  psi -> symmetry    = ORBITA_SYMM_SPHERICAL;
  psi -> coordinates = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_models_plummer_rho;
  psi -> psi = & _orbita_models_plummer_psi;
  psi -> f   = & _orbita_models_plummer_f;
  psi -> Mr  = & _orbita_models_plummer_Mr;
  psi -> Mt  = & _orbita_models_plummer_Mt;

  psi -> is_variable = 0;
  psi -> init_param  = NULL;
  psi -> param_func  = NULL;

  psi -> is_composite = 0;

  psi -> kill = & _orbita_models_plummer_kill;

  return psi;
}

struct orbita_potential *
orbita_potential_plummer_var(double rho,
                             double a)
{
  orbita_err("orbita_potential_homogeneous_sphere_var: Not implemented.");
  return NULL;
}


/* pseudo-isothermal sphere */

double _orbita_models_pseudo_isothermal_rho(const void * par, const double * pos)
{
  double Vc_sq = *((double *)par), Rc_sq = *((double *)par + 1);
  double R_sq = SphRadSq(pos), R = sqrt(R_sq);
  return (Vc_sq / (4. * Pi)) * ((3. * Rc_sq + R_sq) / Sq(Rc_sq + R_sq));
}

double _orbita_models_pseudo_isothermal_psi(const void * par, const double * pos)
{
  double Vc_sq = *((double *)par), Rc_sq = *((double *)par + 1);
  double R_sq = SphRadSq(pos);
  return 0.5 * Vc_sq * log(1. + R_sq / Rc_sq);
}

int _orbita_models_pseudo_isothermal_f(const void * par, const double * pos, double * F)
{
  double Vc_sq = *((double *)par), Rc_sq = *((double *)par + 1);
  double R_sq = SphRadSq(pos), R = sqrt(R_sq);
  double Fr_R = -Vc_sq / Rc_sq / (1. + R_sq / Rc_sq);
  if(! isfinite(Fr_R)) Fr_R = 0;

  F[0] = Fr_R * pos[0],
  F[1] = Fr_R * pos[1],
  F[2] = Fr_R * pos[2];

  return 0;
}

double _orbita_models_pseudo_isothermal_Mr(const void * par, const double R)
{
  double Vc_sq = *((double *)par), Rc_sq = *((double *)par + 1);
  double R_sq = Sq(R);
  return Vc_sq * (R * R_sq / Rc_sq) / (1. + R_sq / Rc_sq);
}

double _orbita_models_pseudo_isothermal_Mt(const void * par)
{
  return 1. / 0.; // Boom!
}

int
_orbita_models_pseudo_isothermal_kill(struct orbita_potential * psi)
{
  if(psi -> is_variable)
    free(psi -> init_param);
  free(psi -> param);
  return 0;
}

struct orbita_potential * orbita_potential_pseudo_isothermal(double Vc, double Rc)
{
  struct orbita_potential * psi =
    (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> param = (void *) malloc(2 * sizeof(double));
  *((double *)(psi -> param))     = Sq(Vc);
  *((double *)(psi -> param) + 1) = Sq(Rc);

  psi -> is_composite = 0;

  psi -> symmetry    = ORBITA_SYMM_SPHERICAL;
  psi -> coordinates = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_models_pseudo_isothermal_rho;
  psi -> psi = & _orbita_models_pseudo_isothermal_psi;
  psi -> f   = & _orbita_models_pseudo_isothermal_f;
  psi -> Mr  = & _orbita_models_pseudo_isothermal_Mr;
  psi -> Mt  = & _orbita_models_pseudo_isothermal_Mt;

  psi -> is_variable = 0;
  psi -> init_param  = NULL;
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_models_pseudo_isothermal_kill;

  return psi;
}
