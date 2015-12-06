
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "orbit.h"
#include "integrator.h"
#include "potential.h"
#include "analyzer.h"

struct orbita_orbit *
orbita_orbit_alloc(const int N_size)
{
  struct orbita_orbit * orb = (struct orbita_orbit *)
    malloc (sizeof (struct orbita_orbit));

  orb -> pts    = (double *) malloc(sizeof(double) * N_size * 6),
  orb -> N_size = N_size,
  orb -> N_pts  = 0;

  return orb;
}

int
orbita_orbit_empty(struct orbita_orbit * orb)
{
  int I_pt, N_len = orb -> N_size * 6;
  for(I_pt = 0; I_pt < N_len; ++ I_pt)
    (orb -> pts)[I_pt] = 0.;
  orb -> N_pts = 0;
  return 0;
}

struct orbita_orbit *
orbita_orbit_integrate(struct orbita_potential *  psi,
                       struct orbita_integrator * intgr,
                       double *                   W_init,
                       double                     T_init,
                       int                        N_steps)
{
  // allocate orbit object
  struct orbita_orbit * orb = (struct orbita_orbit *)
      malloc(sizeof(struct orbita_orbit));
  orb -> pts = (double *) malloc(sizeof(double) * N_steps * 6);

  double * pt = orb -> pts;

  // get ready for integration
  double t = T_init, dt = intgr -> delta_t, wp = intgr -> omega;

  int i_pt;

  // set the initial condition
  for(i_pt = 0; i_pt < 6; ++ i_pt) pt[i_pt] = W_init[i_pt];

  // integrate forwards (and terminate if 'pace' function returns non-zero)
  for(i_pt = 0; i_pt < N_steps - 1; ++ i_pt)
    if ((intgr -> pace)(intgr -> wspace, psi,  pt + 6 * i_pt,
                        dt, t + i_pt * dt, wp, pt + 6 * i_pt + 6)) break;

  // write other information
  orb -> t_init  = t,
  orb -> delta_t = dt;

  orb -> is_rotframe  = intgr -> is_rot_frame,
  orb -> omega        = wp;

  orb -> N_pts  = i_pt,
  orb -> N_size = N_steps;

  return orb;
}

int
orbita_orbit_integrate_noalloc(struct orbita_orbit *      orb,
                               struct orbita_potential *  psi,
                               struct orbita_integrator * intgr,
                               double *                   W_init,
                               double                     T_init,
                               int                        N_steps)
{
  double * pt = orb -> pts;
  int N_size = orb -> N_size;

  if (N_steps > N_size)
    {
      orbita_err ("No enough space for orbit integration.");
      return -1;
    }

  // get ready for integration
  double t = T_init, dt = intgr -> delta_t, wp = intgr -> omega;

  int i_pt;

  // set the initial condition
  for (i_pt = 0; i_pt < 6; ++ i_pt) pt[i_pt] = W_init[i_pt];

  // integrate forwards (and terminate if 'pace' function returns non-zero)
  for(i_pt = 0; i_pt < N_steps - 1; ++ i_pt)
    if ((intgr -> pace) (intgr -> wspace, psi,  pt + 6 * i_pt,
                         dt, t + i_pt * dt, wp, pt + 6 * i_pt + 6))
      break;

  // write other information
  orb -> t_init  = t,
  orb -> delta_t = dt;

  orb -> is_rotframe  = intgr -> is_rot_frame,
  orb -> omega        = wp;

  orb -> N_pts  = i_pt;

  return 0;
}

struct orbita_orbit *
orbita_orbit_integrate_with_analyzer(struct orbita_potential *  psi,
                                     struct orbita_integrator * intgr,
                                     struct orbita_analyzer *   wat,
                                     double *                   W_init,
                                     double                     T_init,
                                     int                        N_steps)
{
  // allocate orbit object
  struct orbita_orbit * orb = (struct orbita_orbit *)
      malloc(sizeof(struct orbita_orbit));
  orb -> pts = (double *) malloc(sizeof(double) * N_steps * 6);
  double * pt = orb -> pts;

  // for integrator
  double t = T_init, dt = intgr -> delta_t, wp = intgr -> omega;

  int i_pt;

  // set the initial condition
  for(i_pt = 0; i_pt < 6; ++ i_pt) pt[i_pt] = W_init[i_pt];

  // integrate forwards, terminate if the analyzer \
     or the integrator return an non-zero value
  for(i_pt = 0; i_pt < N_steps - 1; ++ i_pt)
    if ((intgr -> pace)(intgr -> wspace, psi,  pt + 6 * i_pt,
                        dt, t + i_pt * dt, wp, pt + 6 * i_pt + 6) ||
        orbita_analyzer_peek(wat, pt, & i_pt, N_steps)) break;

  // set other properties
  orb -> t_init  = t,
  orb -> delta_t = dt;

  orb -> is_rotframe  = intgr -> is_rot_frame,
  orb -> omega        = wp;

  orb -> N_pts  = i_pt,
  orb -> N_size = N_steps;

  return orb;
}

int
orbita_orbit_integrate_with_analyzer_noalloc(struct orbita_orbit *      orb,
                                             struct orbita_potential *  psi,
                                             struct orbita_integrator * intgr,
                                             struct orbita_analyzer *   wat,
                                             double *                   W_init,
                                             double                     T_init,
                                             int                        N_steps)
{
  double * pt = orb -> pts;
  int N_size = orb -> N_size;

  if (N_steps > N_size)
    {
      orbita_err ("No enough space for orbit integration.");
      return -1;
    }

  // for integrator
  double t = T_init, dt = intgr -> delta_t, wp = intgr -> omega;

  int i_pt;

  // set the initial condition
  for(i_pt = 0; i_pt < 6; ++ i_pt) pt[i_pt] = W_init[i_pt];

  // integrate forwards, terminate if the analyzer \
     or the integrator return an non-zero value
  for(i_pt = 0; i_pt < N_steps - 1; ++ i_pt)
    if ((intgr -> pace)(intgr -> wspace, psi,  pt + 6 * i_pt,
                        dt, t + i_pt * dt, wp, pt + 6 * i_pt + 6) ||
        orbita_analyzer_peek(wat, pt, & i_pt, N_steps))
      break;

  // set other properties
  orb -> t_init  = t,
  orb -> delta_t = dt;

  orb -> is_rotframe  = intgr -> is_rot_frame,
  orb -> omega        = wp;

  orb -> N_pts  = i_pt;

  return 0;
}

int
orbita_orbit_free(struct orbita_orbit * orb)
{
  free(orb -> pts);
  free(orb);
  return 0;
}

int
orbita_orbit_print(struct orbita_orbit * orb)
{
  int       N_pts = orb -> N_pts;
  double *  pt    = orb -> pts;

  // start printing
  int i_pt;
  for(i_pt = 0; i_pt < N_pts; ++ i_pt)
    printf("*\033[33m% 4u\033[0m  \033[36m% .3f\033[0m  % .3f  % .3f  \033[32m% .3f\033[0m  % .3f  % .3f\n",
           i_pt, pt[i_pt * 6 + 0], pt[i_pt * 6 + 1], pt[i_pt * 6 + 2],
           pt[i_pt * 6 + 3], pt[i_pt * 6 + 4], pt[i_pt * 6 + 5]);
}

int
orbita_orbit_savetxt(struct orbita_orbit * orb, const char * filename)
{
  int       N_pts = orb -> N_pts;
  double *  pt    = orb -> pts;

  FILE * fptr = fopen(filename, "w");

  // start printing
  int i_pt;
  for(i_pt = 0; i_pt < N_pts; ++ i_pt)
    fprintf(fptr, "%f\t%f\t%f\t%f\t%f\t%f\n",
           pt[i_pt * 6 + 0], pt[i_pt * 6 + 1], pt[i_pt * 6 + 2],
           pt[i_pt * 6 + 3], pt[i_pt * 6 + 4], pt[i_pt * 6 + 5]);

  fclose(fptr);

  return 0;
}

int
orbita_orbit_savestdf4(struct orbita_orbit * orb, const char * filename)
{
  int       N_pts = orb -> N_pts;
  double *  pt    = orb -> pts;

  double t0 = orb -> t_init,
         dt = orb -> delta_t;

  FILE * fptr = fopen(filename, "wb");

  // convert to float32
  float * arr = (float *) malloc(sizeof(float) * 7 * N_pts);

  int it, iw;
  for(it = 0; it < N_pts; ++ it)
    {
      for(iw = 0; iw < 6; ++ iw)
        arr[it * 7 + iw] = pt[it * 6 + iw];
      arr[it * 7 + 6] = t0 + it * dt;
    }

  fwrite((void *)arr, sizeof(float), 7 * N_pts, fptr);

  free(arr);
  fclose(fptr);

  return 0;
}

/*
int (* pace)(void * wspace, struct orbita_potential * psi, double * Wi,
             double dt, double t, double omega, double * Wf);
*/

/* allocate a new orbit */

/*
struct orbita_orbit *
orbita_orbit_alloc(int N_points)
{

  struct orbita_orbit * orbit = (struct orbita_orbit *) malloc(sizeof(orbita_orbit));

  orbit -> N_max_pts  = N_points;
  orbit -> N_pts      = 0;
  orbit -> orbit_pts  = (double *) malloc(sizeof(double) * N_points * 6);

  return orbit;
}

int
orbita_orbit_free(struct orbita_orbit * orbit)
{
  free(orbit -> orbit_pts);
  free(orbit);

  return 0;
}

int
orbita_orbit_initialize(const double * pos,
                        const double * vel)
{
  orbit_pts[0] = pos[0], orbit_pts[1] = pos[1], orbit_pts[2] = pos[2];
  orbit_pts[3] = vel[3], orbit_pts[4] = vel[4], orbit_pts[5] = vel[5];

  return 0;
}
*/

/*
  quickdump: For diagnostic purposes only!
  Exported files may not be portable on other computers!
*/

/*
int
orbita_orbit_quickdump(struct orbita_orbit *  Orbit,
                       char *                 Filename)
{
  // open file, failproof
  FILE * dump_file = fopen(Filename, "wb");
  if(dump_file == NULL) return -1;

  fwrite((const void *)(Orbit), sizeof(struct orbita_orbit), 1, dump_file);
  fwrite((const void *)(Orbit -> orbit_pts), sizeof(double), 6 * (Orbit -> N_pts), dump_file);

  fclose(dump_file);

  return 0;
}

struct orbita_orbit *
orbita_orbit_quickload(char * Filename)
{
  FILE * dump_file = fopen(Filename, "rb");
  if(dump_file == NULL) return -1;

  struct orbita_orbit * Orbit = (struct orbita_orbit *) malloc(sizeof(struct orbita_orbit));
  fread((void *)Orbit, sizeof(struct orbita_orbit), 1, dump_file);

  Orbit -> orbit_pts = (double *) malloc(sizeof(double) * 6 * (Orbit -> N_max_pts));
  fread((void *)(Orbit -> orbit_pts), sizeof(double), 6 * (Orbit -> N_pts), dump_file);

  fclose(dump_file);

  return 0;
}

int
orbita_orbit_quicklist(struct orbita_orbit *  Orbit,
                       FILE *                 Stream)
{
  int I_ln; double * pt = Orbit -> orbit_pts,

  for(I_ln = 0; I_ln < (Orbit -> N_pts); ++ I_ln, pt += 6)
    fprintf(Stream, "% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f\n",
                     pt[0],  pt[1],  pt[2],  pt[3],  pt[4],  pt[5] );

  return 0
}
*/
