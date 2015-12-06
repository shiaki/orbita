
// In memory of Claudio Abbado

#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "potential.h"
#include "models.h"

#define Sq(X) ((X)*(X))

unsigned int orbita_potential_azimuthal_average_points = 256;

/* Allocate and free potentials */

struct orbita_potential *
orbita_potential_alloc(enum orbita_potential_type Type,
                       const void *               Init_param,
                       int                        Is_variable,
                       int                        (* Param_func)(void const *,    // p(0), init params
                                                                 const double,    // t,    time
                                                                 void *))         // p(t), output params
{
  struct orbita_potential * psi;

  switch(Type)
    {
      case ORBITA_POTENTIAL_POINTMASS:
        {
          if(Is_variable)
            psi = orbita_potential_pointmass_var(
                      *((double *)(Init_param)), Param_func);
          else psi = orbita_potential_pointmass(*((double *)(Init_param)));
          break;
        }

      case ORBITA_POTENTIAL_HOMOGENEOUS_SPHERE:
        {
          if(Is_variable) orbita_err("homogeneous_sphere is not variable yet.");
          else psi = orbita_potential_homogeneous_sphere(
                   *((double *)Init_param), *((double *)Init_param + 1));
          break;
        }

      /*
      case YOUR_CUSTOM_POTENTIAL:
        {
          // do something
          break;
        }
      */

      default:
        orbita_err("orbita_potential_alloc: Undefined potential type.");
        return NULL;
    }

  psi -> type         = Type;

  psi -> is_composite = 0,
  psi -> N_components = 0,
  psi -> components   = NULL;

  psi -> is_variable  = Is_variable,
  psi -> param_func   = Is_variable?Param_func:NULL;

  return psi;
}

struct orbita_potential *
orbita_potential_alloc_composite(int N_components, ...)
{
  struct orbita_potential * psi =
    (struct orbita_potential *) malloc(sizeof(struct orbita_potential));

  psi -> is_composite = 1;
  psi -> N_components = N_components;
  psi -> components   = (struct orbita_potential **)
           malloc(sizeof(struct orbita_potential *) * N_components);

  // now read parameters
  va_list  arglist;
  va_start(arglist, N_components);

  int I_component;
  for(I_component = 0; I_component < N_components; ++ I_component)
    (psi -> components)[I_component] =
              va_arg(arglist, struct orbita_potential *);

  va_end(arglist);

  return psi;
}

int
orbita_potential_free(struct orbita_potential * psi)
{

  // for composite potential, loop over components
  if(psi -> is_composite)
    {
      int I_component;
      for(I_component = 0; I_component < psi -> N_components; ++ I_component)
        orbita_potential_free((psi -> components)[I_component]);
    }

  // for single-component potential, free parameters, and then free potential
  else
    (psi -> kill)(psi);

  free(psi);

  return 0;
}

/* Force, density, potential, rot-curve, resonance freq, etc. */

int
orbita_potential_evaluate_force(struct orbita_potential * psi, const double * pos, double t, double * F)
{
  if(psi -> is_composite)
    {
      double Ft[3];
      F[0] = 0., F[1] = 0., F[2] = 0.;

      int I_component, It;
      for(I_component = 0; I_component < psi -> N_components; ++ I_component)
        {
          //((psi -> components[I_component]) -> f)(psi -> param, pos, Ft);
          It = orbita_potential_evaluate_force(
                   ((psi -> components)[I_component]), pos, t, Ft);
          F[0] += Ft[0], F[1] += Ft[1], F[2] += Ft[2];
        }
    }

  else
    {
      // update the parameters if the potential is variable
      if(psi -> is_variable) (psi -> param_func)(psi -> init_param, t, psi -> param);
      return (psi -> f)(psi -> param, pos, F);
    }
}

double
orbita_potential_evaluate_density(struct orbita_potential * psi, const double * pos, double t)
{
  if(psi -> is_composite)
    {
      double rho_t = 0.; int I_component;
      for(I_component = 0; I_component < psi -> N_components; ++ I_component)
        rho_t += orbita_potential_evaluate_density(
                     (psi -> components)[I_component], pos, t);
        return rho_t;
    }

  else
    {
      if(psi -> is_variable) (psi -> param_func)(psi -> init_param, t, psi -> param);
      return (psi -> rho)(psi -> param, pos);
    }
}

double
orbita_potential_evaluate_potential(struct orbita_potential * psi, const double * pos, double t)
{
  if(psi -> is_composite)
    {
      double psi_t = 0.; int I_component;
      for(I_component = 0; I_component < psi -> N_components; ++ I_component)
        psi_t += orbita_potential_evaluate_potential(
                     (psi -> components)[I_component], pos, t);
        return psi_t;
    }

  else
    {
      if(psi -> is_variable) (psi -> param_func)(psi -> init_param, t, psi -> param);
      return (psi -> psi)(psi -> param, pos);
    }
}

double
orbita_potential_circular_velocity(struct orbita_potential * psi,
                                   double R,
                                   double t)
{
  // iterate over components, if the potential is composite,
  if(psi -> is_composite)
    {
      double V_cir_sq = 0; int I_component;
      for(I_component = 0; I_component < psi -> N_components; ++ I_component)
        V_cir_sq += Sq(orbita_potential_circular_velocity(psi -> components[I_component], R, t));
      return sqrt(V_cir_sq);
    }
  else
    {
      // if the potential is spherical, or axisymmetric, evaluate the radial force directly
      if((psi -> symmetry) | ORBITA_SYMM_SPHERICAL | ORBITA_SYMM_AXISYMM)
        {
          double pos_t[3] = {R, 0., 0.}, F[3], Fr;
          //(psi -> f)(psi -> param, pos_t, F);
          orbita_potential_evaluate_force(psi, pos_t, t, F);
          return sqrt(fabs(F[0]) * R);
        }
      else
        {
          double pos_t[3], F[3], Fr, V_cir = 0;
          double delta_phi = 2. * Pi / (1. + orbita_potential_azimuthal_average_points);

          int I_pt;
          for(I_pt = 0; I_pt < orbita_potential_azimuthal_average_points; ++ I_pt)
            {
              pos_t[0] = R * cos(I_pt * delta_phi),
              pos_t[1] = R * sin(I_pt * delta_phi),
              pos_t[2] = 0.;

              //(psi -> f)(psi -> param, pos_t, F);
              orbita_potential_evaluate_force(psi, pos_t, t, F);
              Fr = F[0] * cos(I_pt * delta_phi) + F[1] * sin(I_pt * delta_phi);
              V_cir += sqrt(Fr * R);
            }

          V_cir /= orbita_potential_azimuthal_average_points;
          return V_cir;
        }
    }
}

// Vect wrap of orbita_potential_circular_velocity
int
orbita_potential_circular_velocity_vect(struct orbita_potential *  psi,
                                        int                        N_pts,
                                        double *                   R,
                                        double *                   V)
{
  int I_pt;
  for(I_pt = 0; I_pt < N_pts; ++ I_pt)
    V[I_pt] = orbita_potential_circular_velocity(psi, R[I_pt], 0.); /* TODO: the parameter T */

  return 0;
}

double
orbita_potential_circular_frequency(struct orbita_potential *  psi,
                                    double                     Rg)
{
  return orbita_potential_circular_velocity(psi, Rg, 0.) / Rg; /* TODO: the parameter T */
}

int
orbita_potential_circular_frequency_vect(struct orbita_potential *  psi,
                                         int                        N_pts,
                                         double *                   Rg,
                                         double *                   Omega)
{
  int I_pt;
  for(I_pt = 0; I_pt < N_pts; ++ I_pt)
    Omega[I_pt] = orbita_potential_circular_velocity(psi, Rg[I_pt], 0.) / Rg[I_pt]; /* TODO: the parameter T */

  return 0;
}

/*
double
orbita_potential_epicyclic_frequency(struct orbita_potential *  psi,
                                     double                     Rg)
{
}

double
orbita_potential_vertical_frequency(struct orbita_potential *  psi,
                                    double                     Rg)
{
}
*/
