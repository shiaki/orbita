
#include <stdlib.h>
#include <stdio.h>

#include <stddef.h>
#include <float.h>
#include <math.h>

#include <stdlib.h>

#include "potential.h"
#include "integrator.h"
#include "analyzer.h"
#include "orbit.h"

#include "runge_kutta.h"

#include "utils.h"

struct orbita_integrator *
orbita_integrator_init(enum orbita_integrator_type  Type,
                       double                       Timestep,
                       int                          Is_rot_frame,
                       double                       Omega)
{
  struct orbita_integrator * intgr;

  // Inertia frame integrators
  if(!Is_rot_frame)
    switch(Type)
      {
        /*
          case ORBITA_INTEGRATOR_EULER:
          intgr -> step = & _orbita_integrator_euler_pace; break;
        */

        case ORBITA_INTEGRATOR_RUNGE_KUTTA_4:
          intgr = orbita_integrator_runge_kutta_4_init(Timestep); break;

        default:
          orbita_err("Invalid integrator type."); return NULL;
      }

  // rotating frame integrators
  else
    switch(Type)
      {
        case ORBITA_INTEGRATOR_RUNGE_KUTTA_4:
          intgr = orbita_integrator_runge_kutta_4_rot_init(
                    Timestep, Omega); break;

        default:
          orbita_err("Invalid integrator type."); return NULL;
      }

  intgr -> type         = Type,
  intgr -> is_adapt_dt  = 0;

  return intgr;
}

int
orbita_integrator_free(struct orbita_integrator * intgr)
{
  (intgr -> kill)(intgr);
  free(intgr);
  return 0;
}

int
orbita_integrator_optimize_timestep(struct orbita_potential *  psi,
                                    struct orbita_integrator * intgr,
                                    struct orbita_orbit *      orb,
                                    double *                   wi,
                                    double                     eps)
{
  // Reserved for init pos, one/two-step integration
  double * opt_dt_wi = wi;
  double   opt_dt_fs[6], opt_dt_hs1[6], opt_dt_hs2[6];

  // timestep to optimize, and t_init
  double   dt = intgr -> delta_t, t0 = orb -> t_init;

  // desired steps, order of the integrator, and desired accuracy
  double   N_steps = 1000.,//(double)(orb -> N_size),
           k       = (intgr -> order);
           //eps     = intgr -> eps;
           orbita_dbg("get order %f, get steps %f", k, N_steps);

  // integrate one step with normal time step
  (intgr -> pace) (intgr -> wspace, psi,
    opt_dt_wi, 2. * dt, t0, intgr -> omega, opt_dt_fs);

  // integrate two steps with half time step
  (intgr -> pace) (intgr -> wspace, psi,
    opt_dt_wi, dt, t0, intgr -> omega, opt_dt_hs1);
  (intgr -> pace) (intgr -> wspace, psi,
    opt_dt_hs1, dt, t0 + dt, intgr -> omega, opt_dt_hs2);

  // calculate error
  #define Sq(X) ((X) * (X))
  int it; double w0 = 0, dw = 0;
  for(it = 0; it < 6; ++ it)
    w0 += Sq(opt_dt_wi[it]), dw += Sq(opt_dt_hs2[it] - opt_dt_fs[it]);
  w0 = sqrt(w0), dw = sqrt(dw);
  #undef Sq

  // find h_max
  double dt_opt = pow(2. * (pow(2., k) - 1.) * (1. / N_steps)
                      * (eps * w0 / dw), 1. / k) * dt;

  // use optimized dt
  printf("Timestep: optimal %e, current %f\n", dt_opt, dt);
  while(fabs(dt) > dt_opt) dt /= 2.;
  intgr -> delta_t = dt;
  printf("          optimized timestep: %f\n", dt);

  return 0;
}

int
orbita_integrator_update_timestep(struct orbita_integrator * intgr,
                                  double delta_t)
{
  intgr -> delta_t = delta_t;
  return 0;
}
