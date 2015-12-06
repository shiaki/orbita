
/*
  TODO

    1. More generic case, launching towards negative X, or positive X
    2. More generic, searching along X axis or Y
    3. Allow mul2 orbits searching
*/

#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include <math.h>
#include <float.h>

#include "periodic_seeker.h"

#include "utils.h"
#include "potential.h"
#include "analyzer.h"

#define Sq(X) ((X)*(X))
#define GR (0.61803398875)

// Helper function to find vphi
double _pseeker_x_vphi(struct orbita_potential * psi,
                       double E_J, double omega_sq, double R)
{
  double pos_t[3] = {0., R, 0.};
  double vsq_t = 2. * E_J + (R * R) * omega_sq
      - 2. * orbita_potential_evaluate_potential(psi, pos_t, 0.);
  return sqrt(vsq_t);
}

int
_orbita_analyzer_planar_periodic_seeker_EJ_peek(void * wspac,
    double * pt, int * i_pt, const int N_pts)
{
  // unpack the workspace
  struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac *) wspac;

  // first run, make initial condition directly
  if(ws -> I_iter == 0) goto Make_init;

  size_t w0_idx = 6 * (* i_pt),
         w1_idx = 6 * (* i_pt + 1);

  // if no Y-axis crossing, pass
  if(pt[w1_idx] * pt[w0_idx] > 0.) return 0;

  // skip the false start.
  if(* i_pt == 0) return 0;

  if(pt[w1_idx + 3] > 0) return 0; // not positive crossing

  // true Y-axis crossing detected.

  // Find the intersection point using linear interpolation
  double h = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - h) * pt[w0_idx + 1]  + h * pt[w1_idx + 1],
    pxs = (1. - h) * pt[w0_idx + 3]  + h * pt[w1_idx + 3],
    pys = (1. - h) * pt[w0_idx + 4]  + h * pt[w1_idx + 4];

  // Find the phase space distance
  double d_t = Sq((ws -> R_t) - ys) + Sq(pys) + Sq((ws -> V_t) - pxs);

  // delta small enough? exit.
  if(fabs(d_t) < 1e-10) return 1;

  Make_init:;

  double R_t;

  if(ws -> I_iter == 0) // start estimate F_a
    {
    R_t = ws -> R_a;
    orbita_dbg("0th pass, set R_a = %f", R_t);}

  if(ws -> I_iter == 1) // collect F_a, start estimate F_b
    {
    ws -> d_a = d_t, R_t = ws -> R_b;
    orbita_dbg("1st pass, get d_a = %f, set R_b = %f", d_t, R_t);}

  if(ws -> I_iter == 2) // collect F_b, start estimate F_c
    {
    ws -> d_b = d_t,
    R_t = (ws -> R_b) - GR * ((ws -> R_b) - (ws -> R_a)),
    ws -> R_c = R_t;
    orbita_dbg("2nd pass, get d_b = %f, set R_c = %f", d_t, R_t);}

  if(ws -> I_iter == 3) // collect F_c, take left.
    {
    ws -> d_c = d_t,
    R_t = (ws -> R_a) + GR * ((ws -> R_b) - (ws -> R_a)),
    ws -> R_d = R_t, ws -> take_left = 0;
    orbita_dbg("3rd pass, get d_c = %f, set R_d = %f", d_t, R_t);}

  if(ws -> I_iter >= 4)
    {

    if(ws -> take_left) ws -> d_c = d_t,
    orbita_dbg("%u pass, get d_c = %f", ws -> I_iter, d_t);
    else ws -> d_d = d_t,
    orbita_dbg("%u pass, get d_d = %f", ws -> I_iter, d_t);

    orbita_dbg("Ra %f, Rc %f, Rd %f, Rb %f",
        ws -> R_a, ws -> R_c, ws -> R_d, ws -> R_b);
    orbita_dbg("da %f, dc %f, dd %f, db %f",
        ws -> d_a, ws -> d_c, ws -> d_d, ws -> d_b);

  //if(ws -> I_iter >= 4)
    if(ws -> d_c < ws -> d_d)
    {
      ws -> R_b = ws -> R_d, ws -> R_d = ws -> R_c,
      ws -> d_b = ws -> d_d, ws -> d_d = ws -> d_c,
      ws -> R_c = (ws -> R_b) - GR * ((ws -> R_b) - (ws -> R_a)),
      ws -> take_left = 1;
      orbita_dbg("> d_c < d_d, set R_c = %f, take LEFT.", ws -> R_c);
    }
    else
    {
      ws -> R_a = ws -> R_c, ws -> R_c = ws -> R_d,
      ws -> d_a = ws -> d_c, ws -> d_c = ws -> d_d,
      ws -> R_d = (ws -> R_a) + GR * ((ws -> R_b) - (ws -> R_a)),
      ws -> take_left = 0;
      orbita_dbg("> d_c > d_d, set R_d = %f, take RIGHT.", ws -> R_c);
    }

    // TODO> Set R_t Here!!
      }

  // set initial condition
  pt[0] = 0., pt[1] = R_t, pt[2] = 0.,
  pt[3] = -_pseeker_x_vphi(ws -> psi, ws -> E_J, ws -> omega_sq, R_t),
  pt[4] = 0., pt[5] = 0.;

  ws -> R_t = R_t, ws -> V_t = pt[3];

  ws -> I_iter += 1;
  if(ws -> I_iter == ws -> N_iter) return 1;

  * i_pt = -1;

  printf("\n");

  return 0;
}

int
_orbita_analyzer_planar_periodic_seeker_EJ_kill(
    struct orbita_analyzer * wat)
{
  free(wat -> wspac); return 0;
}

struct orbita_analyzer *
orbita_analyzer_planar_periodic_seeker_EJ(struct orbita_potential * psi,
  double omega, double E_J, double R_min, double R_max, int N_iter)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac * ws =
   (struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac *)
   malloc(sizeof(struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac));

  ws -> psi = psi, ws -> omega_sq = omega * omega,
  ws -> I_iter = 0, ws -> N_iter = N_iter,
  ws -> E_J = E_J, ws -> R_a = R_min, ws -> R_b = R_max;
  wat -> wspac = (void *) ws;

  wat -> is_group = 0,
  wat -> N_analyzers = 0,
  wat -> group = NULL;

  wat -> peek = _orbita_analyzer_planar_periodic_seeker_EJ_peek,
  wat -> kill = _orbita_analyzer_planar_periodic_seeker_EJ_kill;

  return wat;
}

/*
  peek function, searching for periodic orbit for a given pos.
*/
int
_orbita_analyzer_planar_periodic_seeker_pos_peek(void * wspac,
    double * pt, int * i_pt, const int N_pts)
{
  // unpack workspace
  struct _orbita_analyzer_planar_periodic_seeker_pos_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_seeker_pos_wspac *) wspac;

  // first run? make initial condition directly,
  if(ws -> I_iter == 0) goto Make_init;

  size_t w0_idx = 6 * (* i_pt),
         w1_idx = 6 * (* i_pt + 1);

  // if no Y-axis crossing, pass
  if(pt[w1_idx] * pt[w0_idx] > 0.) return 0;

  // skip the false start.
  if(* i_pt == 0) return 0;

  // if false crossing, skip
  int is_positive_cross = pt[w1_idx + 3] > 0;
  if((ws -> V_t > 0.) != is_positive_cross) return 0;

  // Now detect true positive crossing.

  // Not this turn? skip.
  ++ (ws -> I_turn);
  if(ws -> I_turn != ws -> N_turn) return 0;

  // Find the intersection point using linear interpolation
  double h = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - h) * pt[w0_idx + 1]  + h * pt[w1_idx + 1],
    pxs = (1. - h) * pt[w0_idx + 3]  + h * pt[w1_idx + 3],
    pys = (1. - h) * pt[w0_idx + 4]  + h * pt[w1_idx + 4];

  // Find the phase space distance
  double d_t = sqrt(Sq((ws -> R) - ys)
             + Sq(pys) + Sq((ws -> V_t) - pxs));

  if(fabs(d_t) < 1e-10) return 1;

  d_t = log(d_t);

  Make_init:;

  orbita_dbg("I_it %u", ws -> I_iter);
  /*
  orbita_dbg("V: a, c, d, b = %f, %f, %f, %f", ws -> V_a,
             ws -> V_c, ws -> V_d, ws -> V_b);
  orbita_dbg("D: a, c, d, b = %f, %f, %f, %f", ws -> d_a,
             ws -> d_c, ws -> d_d, ws -> d_b);
  */

  double V_t;

  if(ws -> I_iter == 0) V_t = ws -> V_a;

  if(ws -> I_iter == 1) V_t = ws -> V_b, ws -> d_a = d_t;

  if(ws -> I_iter == 2)
    V_t = (ws -> V_b) - GR * ((ws -> V_b) - (ws -> V_a)),
    ws -> d_b = d_t, ws -> V_c = V_t;

  if(ws -> I_iter == 3)
    V_t = (ws -> V_a) + GR * ((ws -> V_b) - (ws -> V_a)),
    ws -> d_c = d_t, ws -> V_d = V_t,
    ws -> take_left = 1;

  if(ws -> I_iter >= 4)
    {
      if(ws -> take_left) ws -> d_c = d_t;
      else ws -> d_d = d_t;

  orbita_dbg("Get d = %f", d_t);
  orbita_dbg("V: a, c, d, b = %f, %f, %f, %f", ws -> V_a,
             ws -> V_c, ws -> V_d, ws -> V_b);
  orbita_dbg("D: a, c, d, b = %f, %f, %f, %f", ws -> d_a,
             ws -> d_c, ws -> d_d, ws -> d_b);

      if(ws -> d_c < ws -> d_d)
        {
          ws -> V_b = ws -> V_d, ws -> V_d = ws -> V_c,
          ws -> d_b = ws -> d_d, ws -> d_d = ws -> d_c,
          V_t = (ws -> V_b) - GR * ((ws -> V_b) - (ws -> V_a)),
          ws -> V_c = V_t, ws -> take_left = 1;
        }
      else
        {
          ws -> V_a = ws -> V_c, ws -> V_c = ws -> V_d,
          ws -> d_a = ws -> d_c, ws -> d_c = ws -> d_d,
          V_t = (ws -> V_a) + GR * ((ws -> V_b) - (ws -> V_a)),
          ws -> V_d = V_t, ws -> take_left = 0;
        }
    }

  // set initial condition
  pt[0] = 0., pt[1] = ws -> R, pt[2] = 0.,
  pt[3] = V_t, pt[4] = 0., pt[5] = 0.;

  ws -> V_t = V_t;
  orbita_dbg("set V_t = %f", V_t);

  ws -> I_iter += 1;
  if(ws -> I_iter == ws -> N_iter) return 1;

  * i_pt = -1;

  ws -> I_turn = 0;

  printf("\n");

  return 0;
}

int
_orbita_analyzer_planar_periodic_seeker_pos_kill(
    struct orbita_analyzer * wat)
{
  free(wat -> wspac); return 0;
}

struct orbita_analyzer *
orbita_analyzer_planar_periodic_seeker_pos(double R,
  double V_a, double V_b, int pos_vel, int N_iter, int N_turn)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_planar_periodic_seeker_pos_wspac * ws =
  (struct _orbita_analyzer_planar_periodic_seeker_pos_wspac *)
  malloc(sizeof(struct _orbita_analyzer_planar_periodic_seeker_pos_wspac));

  ws -> I_turn = 0, ws -> N_turn = N_turn,
  ws -> I_iter = 0, ws -> N_iter = N_iter,
  ws -> R = R, ws -> V_a = V_a, ws -> V_b = V_b;
  wat -> wspac = (void *) ws;

  wat -> is_group = 0,
  wat -> N_analyzers = 0,
  wat -> group = NULL;

  wat -> peek = _orbita_analyzer_planar_periodic_seeker_pos_peek,
  wat -> kill = _orbita_analyzer_planar_periodic_seeker_pos_kill;

  return wat;
}


/* quick fox jumps over a lazy dog */

int
_orbita_analyzer_planar_periodic_ranger_pos_peek(void * wspac,
    double * pt, int * i_pt, const int N_pts)
{
  // unpack the workspace
  struct _orbita_analyzer_planar_periodic_ranger_pos_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_pos_wspac *) wspac;

  // first run, go make initial condition directly
  if((ws -> is_init) ++ == 0) goto Make_init;

  // this point and (already integrated) next point
  size_t w0_idx = 6 * (* i_pt),
         w1_idx = 6 * (* i_pt + 1);

  // if no Y-axis crossing, pass
  if(pt[w1_idx] * pt[w0_idx] > 0.) return 0;

  // skip the false start.
  if(* i_pt == 0) return 0;

  // if false crossing, skip
  int is_positive_cross = pt[w1_idx + 3] > 0;
  if(((ws -> V)[ws -> I_pt] > 0.) != is_positive_cross) return 0;

  // Not this turn? skip.
  ++ (ws -> I_turn);
  if(ws -> I_turn != ws -> N_turn) return 0;

  // find the intersection point
  double h = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - h) * pt[w0_idx + 1]  + h * pt[w1_idx + 1],
    pxs = (1. - h) * pt[w0_idx + 3]  + h * pt[w1_idx + 3],
    pys = (1. - h) * pt[w0_idx + 4]  + h * pt[w1_idx + 4];

  // Find the phase space distance
  double d_t = sqrt(Sq((ws -> R) - ys)
             + Sq(pys) + Sq((ws -> V)[ws -> I_pt] - pxs));
  (ws -> D)[ws -> I_pt] = d_t;

  // DBG
  //printf("%f %f\n", (ws -> V)[ws -> I_pt], d_t);

  // go to the next point
  ++ (ws -> I_pt);
  if(ws -> I_pt == ws -> N_pts) return 1;

  Make_init:;

  double V_t = (ws -> V)[ws -> I_pt];

  // set initial condition
  pt[0] = 0., pt[1] = ws -> R, pt[2] = 0.,
  pt[3] = V_t, pt[4] = 0., pt[5] = 0.;

  ws -> I_turn = 0, * i_pt = -1;

  return 0;
}

int
_orbita_analyzer_planar_periodic_ranger_pos_kill(
    struct orbita_analyzer * wat)
{
  struct _orbita_analyzer_planar_periodic_ranger_pos_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_pos_wspac *)
    (wat -> wspac);

  free(ws -> V), free(ws -> D);
  free(wat -> wspac);

  return 0;
}

struct orbita_analyzer *
orbita_analyzer_planar_periodic_ranger_pos(double R, double V_a,
    double V_b, int N_pts, int N_turn)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_planar_periodic_ranger_pos_wspac * ws =
  (struct _orbita_analyzer_planar_periodic_ranger_pos_wspac *)
  malloc(sizeof(struct _orbita_analyzer_planar_periodic_ranger_pos_wspac));
  wat -> wspac = (void *) ws;

  ws -> I_turn = 0, ws -> N_turn = N_turn,
  ws -> N_pts = N_pts, ws -> I_pt = 0,
  ws -> R = R;

  ws -> V = (double *) malloc(sizeof(double) * N_pts),
  ws -> D = (double *) malloc(sizeof(double) * N_pts);

  ws -> is_init = 0;

  int i_pt;
  for(i_pt = 0; i_pt < N_pts; ++ i_pt)
    (ws -> V)[i_pt] = V_a + i_pt * (V_b - V_a) / (N_pts - 1),
    (ws -> D)[i_pt] = 0.;

  wat -> is_group = 0,
  wat -> N_analyzers = 0,
  wat -> group = NULL;

  wat -> peek = _orbita_analyzer_planar_periodic_ranger_pos_peek,
  wat -> kill = _orbita_analyzer_planar_periodic_ranger_pos_kill;

  return wat;
}

int
orbita_analyzer_planar_periodic_ranger_pos_get(
    struct orbita_analyzer * wat, double * V, double * D)
{
  struct _orbita_analyzer_planar_periodic_ranger_pos_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_pos_wspac *)
    (wat -> wspac);

  // copy results
  memcpy(V, ws -> V, sizeof(double) * (ws -> N_pts)),
  memcpy(D, ws -> D, sizeof(double) * (ws -> N_pts));

  return ws -> N_pts;
}

/* Meow, meow, meow, meow, meow... */

int
_orbita_analyzer_planar_periodic_ranger_EJ_peek(void * wspac,
    double * pt, int * i_pt, const int N_pts)
{
  // unpack the workspace
  struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac *) wspac;

  // first run, go make initial condition directly
  if((ws -> is_init) ++ == 0) goto Make_init;

  // this point and (already integrated) next point
  size_t w0_idx = 6 * (* i_pt),
         w1_idx = 6 * (* i_pt + 1);

  //printf("%f %f %f %f %f %f\n", pt[w0_idx], pt[w0_idx + 1],
  //  pt[w0_idx + 2], pt[w0_idx + 3], pt[w0_idx + 4], pt[w0_idx + 5]);

  //orbita_dbg("at A: I_pt = %u", ws -> I_pt);

  // if fly away, re-make initial condition
  double rsq_t = Sq(pt[w0_idx]) + Sq(pt[w0_idx + 1]) + Sq(pt[w0_idx + 2]);
  if(rsq_t > (ws -> r_max_sq) || isnan(rsq_t))
    {(ws -> D)[(ws -> I_pt) ++] = 1. / 0.; goto Make_init;}

  // if no Y-axis crossing, pass
  if(pt[w1_idx] * pt[w0_idx] > 0.) return 0;

  // skip the false start.
  if(* i_pt == 0) return 0;

  // if false crossing, skip
  int is_positive_cross = pt[w1_idx + 3] > 0;
  if((ws -> take_positive) != is_positive_cross) return 0;

  // Not this turn? skip.
  ++ (ws -> I_turn);
  if(ws -> I_turn != ws -> N_turn) return 0;

  // find the intersection point
  double h = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - h) * pt[w0_idx + 1]  + h * pt[w1_idx + 1],
    pxs = (1. - h) * pt[w0_idx + 3]  + h * pt[w1_idx + 3],
    pys = (1. - h) * pt[w0_idx + 4]  + h * pt[w1_idx + 4];

  //orbita_dbg("at B: I_pt = %u", ws -> I_pt);

  // Find the phase space distance
  double d_t = sqrt(Sq((ws -> R)[ws -> I_pt] - ys)
             + Sq(pys) + Sq((ws -> V_t) - pxs));
  (ws -> D)[ws -> I_pt] = d_t;

  //orbita_dbg("at C: I_pt = %u", ws -> I_pt);

  // go to the next point
  ++ (ws -> I_pt);
  if(ws -> I_pt >= ws -> N_pts) return 1;

  //orbita_dbg("at D: I_pt = %u", ws -> I_pt);

  Make_init:;

  double V_t;

  for(;;)
    {
      if((ws -> I_pt) >= (ws -> N_pts)) return 1;
  //orbita_dbg("at E: I_pt = %u", ws -> I_pt);
      V_t = _pseeker_x_vphi(ws -> psi, ws -> E_J,
            ws -> omega_sq, (ws -> R)[ws -> I_pt]);
  //orbita_dbg("at F: I_pt = %u", ws -> I_pt);
      if(!isnan(V_t)) break;
      else (ws -> D)[(ws -> I_pt) ++] = sqrt(-1.);
  //orbita_dbg("at G: I_pt = %u", ws -> I_pt);
      //orbita_dbg("Energy not allowed");
      //if((ws -> I_pt) >= (ws -> N_pts)) return 1;
  //orbita_dbg("at H: I_pt = %u", ws -> I_pt);
    }

  //orbita_dbg("at I: I_pt = %u", ws -> I_pt);
  // set initial condition
  pt[0] = 0., pt[1] = (ws -> R)[ws -> I_pt], pt[2] = 0.,
  pt[3] = V_t, pt[4] = 0., pt[5] = 0.;

  ws -> V_t = V_t;

  ws -> I_turn = 0, * i_pt = -1;

  //orbita_dbg("V_t = %f", V_t);

  return 0;
}

int
_orbita_analyzer_planar_periodic_ranger_EJ_kill(
    struct orbita_analyzer * wat)
{
  struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac *) (wat -> wspac);

  free(ws -> R), free(ws -> D);
  free(wat -> wspac);

  return 0;
}

struct orbita_analyzer *
orbita_analyzer_planar_periodic_ranger_EJ(struct orbita_potential * psi,
    double omega, double E_J, double R_a, double R_b, double r_max,
    int N_pts, int N_turn)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac * ws =
  (struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac *)
  malloc(sizeof(struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac));
  wat -> wspac = (void *) ws;

  ws -> I_turn = 0, ws -> N_turn = N_turn,
  ws -> N_pts = N_pts, ws -> I_pt = 0,
  ws -> E_J = E_J, ws -> omega_sq = Sq(omega),
  ws -> take_positive = 1, ws -> is_init = 0,
  ws -> r_max_sq = Sq(r_max);

  ws -> psi = psi;

  ws -> R = (double *) malloc(sizeof(double) * N_pts),
  ws -> D = (double *) malloc(sizeof(double) * N_pts);

  int i_pt;
  for(i_pt = 0; i_pt < N_pts; ++ i_pt)
    (ws -> R)[i_pt] = R_a + i_pt * (R_b - R_a) / (N_pts - 1),
    (ws -> D)[i_pt] = 0.;

  wat -> is_group = 0,
  wat -> N_analyzers = 0,
  wat -> group = NULL;

  wat -> peek = _orbita_analyzer_planar_periodic_ranger_EJ_peek,
  wat -> kill = _orbita_analyzer_planar_periodic_ranger_EJ_kill;

  return wat;
}

int
orbita_analyzer_planar_periodic_ranger_EJ_get(
    struct orbita_analyzer * wat, double * R, double * D)
{
  struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac * ws =
    (struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac *)
    (wat -> wspac);

  //orbita_dbg("Called.");

  // copy results
  memcpy(R, ws -> R, sizeof(double) * (ws -> N_pts)),
  memcpy(D, ws -> D, sizeof(double) * (ws -> N_pts));

  return ws -> N_pts;
}

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = */

int
_orbita_analyzer_stopwatch_peek(void * wspac,
    double * pt, int * i_pt, const int N_pts)
{
  // unpack the workspace
  struct _orbita_analyzer_stopwatch_wspac * ws =
    (struct _orbita_analyzer_stopwatch_wspac *) wspac;

  // this point and (already integrated) next point
  size_t w0_idx = 6 * (* i_pt), w1_idx = 6 * (* i_pt + 1);

  // if fly away, terminate and return bad PSD.
  double rsq_t = Sq(pt[w0_idx]) + Sq(pt[w0_idx + 1]) + Sq(pt[w0_idx + 2]);
  if(rsq_t > (ws -> r_max_sq) || isnan(rsq_t)) {(ws -> D) = 1. / 0.; return 1;}

  // if no Y-axis crossing, pass
  if(pt[w1_idx] * pt[w0_idx] > 0.) return 0;

  // skip the start.
  if(* i_pt == 0)
    {
      int i_w;
      for (i_w = 0; i_w < 6; ++ i_w)
        (ws -> W_i)[i_w] = pt[w0_idx + i_w],
        (ws -> W_f)[i_w] = 0. / 0.;
      ws -> take_positive = pt[w0_idx + 3] >= 0.;
      return 0;
    }

  // if false crossing, skip
  int is_positive_cross = pt[w1_idx + 3] >= 0.;
  if(!(ws -> take_positive) != !is_positive_cross) return 0;

  // Not this turn? skip.
  ++ (ws -> I_turn);
  if(ws -> I_turn != ws -> N_turn) return 0;

  // find the intersection point
  double h = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - h) * pt[w0_idx + 1]  + h * pt[w1_idx + 1],
    zs  = (1. - h) * pt[w0_idx + 2]  + h * pt[w1_idx + 2],
    pxs = (1. - h) * pt[w0_idx + 3]  + h * pt[w1_idx + 3],
    pys = (1. - h) * pt[w0_idx + 4]  + h * pt[w1_idx + 4],
    pzs = (1. - h) * pt[w0_idx + 5]  + h * pt[w1_idx + 5];

  // Find the phase space distance
  double d_t = sqrt(Sq((ws -> W_i)[1] - ys) + Sq((ws -> W_i)[2] - zs)
                  + Sq((ws -> W_i)[3] - pxs) + Sq((ws -> W_i)[4] - pys)
                  + Sq((ws -> W_i)[5] - pzs));
  (ws -> D) = d_t;

  (ws -> W_f)[0] = 0.,  (ws -> W_f)[1] = ys,  (ws -> W_f)[2] = zs,
  (ws -> W_f)[3] = pxs, (ws -> W_f)[4] = pys, (ws -> W_f)[5] = pzs;

  // reset stopwatch.
  ws -> I_turn = 0;

  // stop.
  return 1;
}

int
_orbita_analyzer_stopwatch_kill(struct orbita_analyzer * wat)
{
  free(wat -> wspac);
  return 0;
}

struct orbita_analyzer *
orbita_analyzer_stopwatch(const int N_turn, const double R_max)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_stopwatch_wspac * ws =
    (struct _orbita_analyzer_stopwatch_wspac *)
      malloc(sizeof(struct _orbita_analyzer_stopwatch_wspac));
  wat -> wspac = (void *) ws;

  ws -> I_turn        = 0,
  ws -> N_turn        = N_turn,
  //ws -> take_positive = 1,
  ws -> r_max_sq      = Sq(R_max),
  ws -> D             = 1. / 0.;

  wat -> is_group     = 0,
  wat -> N_analyzers  = 0,
  wat -> group        = NULL;

  wat -> peek = _orbita_analyzer_stopwatch_peek,
  wat -> kill = _orbita_analyzer_stopwatch_kill;

  return wat;
}

double
orbita_analyzer_stopwatch_get_D(struct orbita_analyzer * wat)
{
  // unpack the workspace
  struct _orbita_analyzer_stopwatch_wspac * ws =
    (struct _orbita_analyzer_stopwatch_wspac *) (wat -> wspac);

  double rval = ws -> D;
  //ws -> D = DBL_MAX;

  return rval;
}

int
orbita_analyzer_stopwatch_get_Wi(struct orbita_analyzer * wat, double * Wi)
{
  // unpack the workspace
  struct _orbita_analyzer_stopwatch_wspac * ws =
    (struct _orbita_analyzer_stopwatch_wspac *) (wat -> wspac);

  Wi[0] = (ws -> W_i)[0], Wi[1] = (ws -> W_i)[1], Wi[2] = (ws -> W_i)[2],
  Wi[3] = (ws -> W_i)[3], Wi[4] = (ws -> W_i)[4], Wi[5] = (ws -> W_i)[5];

  return 0;
}

int
orbita_analyzer_stopwatch_get_Wf(struct orbita_analyzer * wat, double * Wf)
{
  // unpack the workspace
  struct _orbita_analyzer_stopwatch_wspac * ws =
    (struct _orbita_analyzer_stopwatch_wspac *) (wat -> wspac);

  Wf[0] = (ws -> W_f)[0], Wf[1] = (ws -> W_f)[1], Wf[2] = (ws -> W_f)[2],
  Wf[3] = (ws -> W_f)[3], Wf[4] = (ws -> W_f)[4], Wf[5] = (ws -> W_f)[5];

  return 0;
}
