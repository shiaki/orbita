
#ifndef PERIODIC_SEEKER_H

//#include "analyzer.h"
//#include "potential.h"

struct orbita_potential;
struct orbita_analyzer;

struct _orbita_analyzer_planar_periodic_seeker_EJ_wspac
{
  struct orbita_potential * psi;  // assiciated potential
  double omega_sq;                // pattern speed squared

  int I_iter, N_iter;

  double E_J;
  double R_t, V_t;

  int take_left;

  double R_a, R_b, R_c, R_d;
  double d_a, d_b, d_c, d_d;
};

struct _orbita_analyzer_planar_periodic_seeker_pos_wspac
{
  int take_left;

  int I_turn, N_turn;
  int I_iter, N_iter;

  double R, V_t;

  // for golden method.
  double V_a, V_c, V_d, V_b;
  double d_a, d_c, d_d, d_b;
};

struct _orbita_analyzer_planar_periodic_ranger_pos_wspac
{
  int N_pts;
  int I_pt;

  int N_turn, I_turn;
  int is_init;

  double R;
  double * V, * D;
};

struct _orbita_analyzer_planar_periodic_ranger_EJ_wspac
{
  int N_pts;
  int I_pt;

  struct orbita_potential * psi;

  int take_positive;

  double r_max_sq;

  int N_turn, I_turn;
  int is_init;

  double E_J, omega_sq;
  double * R, * D;

  double V_t;
};

struct _orbita_analyzer_stopwatch_wspac
{
  int take_positive;

  double r_max_sq;

  double W_i[6], W_f[6], D;

  int N_turn, I_turn;
};


struct orbita_analyzer *
orbita_analyzer_planar_periodic_seeker_EJ(struct orbita_potential * psi,
  double omega, double E_J, double R_min, double R_max, int N_iter);

struct orbita_analyzer *
orbita_analyzer_planar_periodic_seeker_pos(double R,
  double V_a, double V_b, int pos_vel, int N_iter, int N_turn);

struct orbita_analyzer *
orbita_analyzer_planar_periodic_ranger_EJ(struct orbita_potential * psi,
    double omega, double E_J, double R_a, double R_b, double r_max,
    int N_pts, int N_turn);

struct orbita_analyzer *
orbita_analyzer_planar_periodic_ranger_pos(double R, double V_a,
    double V_b, int N_pts, int N_turn);

int orbita_analyzer_planar_periodic_ranger_pos_get(
    struct orbita_analyzer * wat, double * V, double * D);

int orbita_analyzer_planar_periodic_ranger_EJ_get(
    struct orbita_analyzer * wat, double * R, double * D);

double orbita_analyzer_stopwatch_get_D(struct orbita_analyzer * wat);

struct orbita_analyzer *
orbita_analyzer_stopwatch(const int N_turn, const double R_max);

int orbita_analyzer_stopwatch_get_Wi(struct orbita_analyzer * wat, double * Wi);
int orbita_analyzer_stopwatch_get_Wf(struct orbita_analyzer * wat, double * Wf);

#define PERIODIC_SEEKER_H
#endif
