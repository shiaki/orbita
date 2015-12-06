
#ifndef ORBITA_ORBIT_H

struct orbita_potential;
struct orbita_integrator;
struct orbita_orbit;
struct orbita_analyzer;

struct orbita_orbit
{
  int N_pts, N_size;

  /*
  enum orbita_coordinate     coordinates;
  struct orbita_potential *  potential;
  struct orbita_integrator * integrator;
  */

  double    t_init;
  double    delta_t;

  int       is_rotframe;
  double    omega;

  //double    init_phase;

  double *  pts;
};

int orbita_orbit_empty(struct orbita_orbit *);

struct orbita_orbit *
orbita_orbit_alloc(const int N_size);

struct orbita_orbit *
orbita_orbit_integrate(struct orbita_potential * psi,
  struct orbita_integrator * intgr, double * W_init,
  double T_init, int N_steps);

int
orbita_orbit_integrate_noalloc(struct orbita_orbit * orb,
    struct orbita_potential * psi,
    struct orbita_integrator * intgr, double * W_init,
    double T_init, int N_steps);

struct orbita_orbit *
orbita_orbit_integrate_with_analyzer(struct orbita_potential * psi,
  struct orbita_integrator * intgr, struct orbita_analyzer * wat,
  double * W_init, double T_init, int N_steps);

int
orbita_orbit_integrate_with_analyzer_noalloc(struct orbita_orbit * orb,
  struct orbita_potential * psi, struct orbita_integrator * intgr,
  struct orbita_analyzer * wat,
  double * W_init, double T_init, int N_steps);

int orbita_orbit_free(struct orbita_orbit *);

int orbita_orbit_print(struct orbita_orbit *);
int orbita_orbit_savetxt(struct orbita_orbit *, const char *);
int orbita_orbit_savestdf4(struct orbita_orbit *, const char *);

#define ORBITA_ORBIT_H
#endif
