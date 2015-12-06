
#ifndef ORBITA_INTEGRATOR_H

struct orbita_potential;
struct orbita_integrator;
struct orbita_orbit;

enum orbita_integrator_type
{

  ORBITA_INTEGRATOR_EULER,

  ORBITA_INTEGRATOR_RUNGE_KUTTA_4,
  ORBITA_INTEGRATOR_RUNGE_KUTTA_8,

  ORBITA_INTEGRATOR_LEAPFROG
};

struct orbita_integrator
{
  /* basic information */
  enum orbita_integrator_type type;
  int  order;

  /* time step */
  double  delta_t;       // Basic step size for integration

  int     is_adapt_dt;   // if we should use adaptive timestep
  double  eps;           // desired accuracy for adaptive timestep

  /* frame rotation */
  int     is_rot_frame;
  double  omega;

  int     is_initialized;
  void *  wspace;

  int (* pace)(void * wspace, struct orbita_potential * psi, double * Wi,
               double dt, double t, double omega, double * Wf);

  int (* kill)(struct orbita_integrator * intgr);

};

struct orbita_integrator *
orbita_integrator_init(enum orbita_integrator_type  Type, double Timestep,
    int Is_rot_frame, double Omega);

int orbita_integrator_free(struct orbita_integrator * intgr);

int orbita_integrator_optimize_timestep(struct orbita_potential *,
    struct orbita_integrator *, struct orbita_orbit *, double *, double);

#define ORBITA_INTEGRATOR_H
#endif
