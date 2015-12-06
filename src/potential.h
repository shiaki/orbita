
#ifndef ORBITA_POTENTIAL_H

enum orbita_potential_type
{
  ORBITA_POTENTIAL_NULL,
  ORBITA_POTENTIAL_POINTMASS,
  ORBITA_POTENTIAL_HOMOGENEOUS_SPHERE
};

enum orbita_coordinate
{
  ORBITA_COORD_CARTESIAN    = 0x00,
  ORBITA_COORD_CYLINDRICAL  = 0x01,
  ORBITA_COORD_SPHERICAL    = 0x02
};

enum orbita_symmetry
{
  ORBITA_SYMM_NONE          = 0x00,

  ORBITA_SYMM_PLANAR        = 0x01,
  ORBITA_SYMM_TRIPLANAR     = 0x02,

  ORBITA_SYMM_TWOFOLD       = 0x04,
  ORBITA_SYMM_AXISYMM       = 0x08,

  ORBITA_SYMM_SPHERICAL     = 0x10
};

#define IS_AXISYMM(psi) ((((psi) -> symmetry & 0x08)?1:0))

struct orbita_potential
{
  // For composite potentials
  int    is_composite;
  int    N_components;

  struct orbita_potential ** components;

  // For single-component potentials
  enum orbita_potential_type  type;
  enum orbita_symmetry        symmetry;
  enum orbita_coordinate      coordinates;

  // Potential parameters
  void * param;

  // For variable parameters
  int     is_variable;
  void *  init_param;
  int  (* param_func)(void const *,   // p(0), init param
                      double,         // t   , time
                      void *);        // p(t), param at time t

  // Func ptrs for potential, force and density
  double (* rho)(void const *, double const *);
  double (* psi)(void const *, double const *);
  int    (* f  )(void const *, double const *, double *);

  // Enclosed mass and total mass
  double (* Mr )(const void *, double);
  double (* Mt )(const void *);

  // Frequencies
  double (* p2R)(const void *, double);
  double (* p2Z)(const void *, double);

  int (* kill)(struct orbita_potential *);
};

struct orbita_potential * orbita_potential_alloc(enum orbita_potential_type,
    const void *, int, int (*)(void const *, const double, void *));
struct orbita_potential * orbita_potential_alloc_composite(int, ...);

int orbita_potential_free(struct orbita_potential *);

int orbita_potential_evaluate_force(struct orbita_potential * psi, const double * pos, double t, double * F);
double orbita_potential_evaluate_density(struct orbita_potential * psi, const double * pos, double t);
double orbita_potential_evaluate_potential(struct orbita_potential * psi, const double * pos, double t);
double orbita_potential_circular_velocity(struct orbita_potential * psi, double R, double t);
int orbita_potential_circular_velocity_vect(struct orbita_potential *  psi, int N_pts, double * R, double * V);
double orbita_potential_circular_frequency(struct orbita_potential *  psi, double Rg);
int orbita_potential_circular_frequency_vect(struct orbita_potential *  psi, int N_pts, double * Rg, double * Omega);

#define ORBITA_POTENTIAL_H
#endif
