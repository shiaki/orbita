
#include "stdio.h"
#include "stdlib.h"

#include "utils.h"
#include "potential.h"
#include "integrator.h"
#include "runge_kutta.h"
/* 4th order Runge-Kutta method */

double  rk4_k1[6], rk4_k2[6], rk4_k3[6], rk4_k4[6],
        rk4_p2[6], rk4_p3[6], rk4_p4[6];

// integrate one step ahead
int
orbita_integrator_runge_kutta_4_lite_step(//struct orbita_integrator * intgr,
                                          struct orbita_potential *  psi,
                                          double *                   wi,
                                          double                     dt,
                                          double                     t,
                                          double                     omega,
                                          double *                   wf)
{
  int stat, it;

  // k1 = f(ti, wi)
  stat = orbita_potential_evaluate_force(psi, wi, t, rk4_k1 + 3);
  for(it = 0; it < 3; ++ it) rk4_k1[it] = wi[it + 3];

  // k2 = f(ti + 0.5 dt, wi + 0.5 * dt * k1)
  for(it = 0; it < 6; ++ it) rk4_p2[it] = wi[it] + 0.5 * dt * rk4_k1[it];
  stat = orbita_potential_evaluate_force(psi, rk4_p2, t + 0.5 * dt, rk4_k2 + 3);
  for(it = 0; it < 3; ++ it) rk4_k2[it] = rk4_p2[it + 3];

  // k3 = f(ti + 0.5 dt, wi + 0.5 * dt * k2)
  for(it = 0; it < 6; ++ it) rk4_p3[it] = wi[it] + 0.5 * dt * rk4_k2[it];
  stat = orbita_potential_evaluate_force(psi, rk4_p3, t + 0.5 * dt, rk4_k3 + 3);
  for(it = 0; it < 3; ++ it) rk4_k3[it] = rk4_p3[it + 3];

  // k4 = f(ti + dt, wi + dt * h3)
  for(it = 0; it < 6; ++ it) rk4_p4[it] = wi[it] + dt * rk4_k3[it];
  stat = orbita_potential_evaluate_force(psi, rk4_p4, t + dt, rk4_k4 + 3);
  for(it = 0; it < 3; ++ it) rk4_k4[it] = rk4_p4[it + 3];

  // wf = wi + (1/6)(k1 + 2 k2 + 2 k3 + k4)
  for(it = 0; it < 6; ++ it)
    wf[it] = wi[it] + (dt / 6.) * (rk4_k1[it]
                + 2. * rk4_k2[it] + 2. * rk4_k3[it] + rk4_k4[it]);

  return 0;
}

int
orbita_integrator_runge_kutta_4_step(void *                     wspace,
                                     struct orbita_potential *  psi,
                                     double *                   wi,
                                     double                     dt,
                                     double                     t,
                                     double                     omega,
                                     double *                   wf)
{
  // unpack workspace
  struct orbita_integrator_runge_kutta_4_workspace * ws =
      (struct orbita_integrator_runge_kutta_4_workspace *)wspace;
  double * k1 = ws -> K1, * k2 = ws -> K2, * k3 = ws -> K3, * k4 = ws -> K4;
  double * p2 = ws -> P2, * p3 = ws -> P3, * p4 = ws -> P4;

  int stat, it;

  // k1 = f(ti, wi)
  stat = orbita_potential_evaluate_force(psi, wi, t, k1 + 3);
  for(it = 0; it < 3; ++ it) k1[it] = wi[it + 3];

  // k2 = f(ti + 0.5 dt, wi + 0.5 * dt * k1)
  for(it = 0; it < 6; ++ it) p2[it] = wi[it] + 0.5 * dt * k1[it];
  stat = orbita_potential_evaluate_force(psi, p2, t + 0.5 * dt, k2 + 3);
  for(it = 0; it < 3; ++ it) k2[it] = p2[it + 3];

  // k3 = f(ti + 0.5 dt, wi + 0.5 * dt * k2)
  for(it = 0; it < 6; ++ it) p3[it] = wi[it] + 0.5 * dt * k2[it];
  stat = orbita_potential_evaluate_force(psi, p3, t + 0.5 * dt, k3 + 3);
  for(it = 0; it < 3; ++ it) k3[it] = p3[it + 3];

  // k4 = f(ti + dt, wi + dt * h3)
  for(it = 0; it < 6; ++ it) p4[it] = wi[it] + dt * k3[it];
  stat = orbita_potential_evaluate_force(psi, p4, t + dt, k4 + 3);
  for(it = 0; it < 3; ++ it) k4[it] = p4[it + 3];

  // wf = wi + (1/6)(k1 + 2 k2 + 2 k3 + k4)
  for(it = 0; it < 6; ++ it)
    wf[it] = wi[it] + (dt / 6.) * (k1[it] + 2. * k2[it]
                 + 2. * k3[it] + k4[it]);

  return 0;
}

int
orbita_integrator_runge_kutta_4_rot_step(void *                     wspace,
                                         struct orbita_potential *  psi,
                                         double *                   wi,
                                         double                     dt,
                                         double                     t,
                                         double                     omega,
                                         double *                   wf)
{
  // unpack workspace
  struct orbita_integrator_runge_kutta_4_workspace * ws =
      (struct orbita_integrator_runge_kutta_4_workspace *)wspace;
  double * k1 = ws -> K1, * k2 = ws -> K2, * k3 = ws -> K3, * k4 = ws -> K4;
  double * p2 = ws -> P2, * p3 = ws -> P3, * p4 = ws -> P4;

  int    stat, it;
  double wp_sq = omega * omega;

  // k1 = f(ti, wi)
  stat = orbita_potential_evaluate_force(psi, wi, t, k1 + 3);
  k1[3] += 2. * omega * wi[4] + wp_sq * wi[0],
  k1[4] -= 2. * omega * wi[3] - wp_sq * wi[1];
  for(it = 0; it < 3; ++ it) k1[it] = wi[it + 3];

  //orbita_dbg("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f", k1[0], k1[1], k1[2], k1[3], k1[4], k1[5]);

  // k2 = f(ti + 0.5 dt, wi + 0.5 * dt * k1)
  for(it = 0; it < 6; ++ it) p2[it] = wi[it] + 0.5 * dt * k1[it];
  stat = orbita_potential_evaluate_force(psi, p2, t + 0.5 * dt, k2 + 3);
  k2[3] += 2. * omega * p2[4] + wp_sq * p2[0],
  k2[4] -= 2. * omega * p2[3] - wp_sq * p2[1];
  for(it = 0; it < 3; ++ it) k2[it] = p2[it + 3];

  // k3 = f(ti + 0.5 dt, wi + 0.5 * dt * k2)
  for(it = 0; it < 6; ++ it) p3[it] = wi[it] + 0.5 * dt * k2[it];
  stat = orbita_potential_evaluate_force(psi, p3, t + 0.5 * dt, k3 + 3);
  k3[3] += 2. * omega * p3[4] + wp_sq * p3[0],
  k3[4] -= 2. * omega * p3[3] - wp_sq * p3[1];
  for(it = 0; it < 3; ++ it) k3[it] = p3[it + 3];

  // k4 = f(ti + dt, wi + dt * h3)
  for(it = 0; it < 6; ++ it) p4[it] = wi[it] + dt * k3[it];
  stat = orbita_potential_evaluate_force(psi, p4, t + dt, k4 + 3);
  k4[3] += 2. * omega * p4[4] + wp_sq * p4[0],
  k4[4] -= 2. * omega * p4[3] - wp_sq * p4[1];
  for(it = 0; it < 3; ++ it) k4[it] = p4[it + 3];

  // wf = wi + (1/6)(k1 + 2 k2 + 2 k3 + k4)
  for(it = 0; it < 6; ++ it)
    wf[it] = wi[it] + (dt / 6.) * (k1[it] + 2. * k2[it]
                 + 2. * k3[it] + k4[it]);

  return 0;
}

int
orbita_integrator_runge_kutta_4_kill(struct orbita_integrator * intgr)
{
  if(intgr -> wspace != NULL) free(intgr -> wspace);
  else orbita_err("Cannot free this integrator. Already removed?");
  return 0;
}

struct orbita_integrator *
orbita_integrator_runge_kutta_4_init(double dt)
{
  struct orbita_integrator * intgr =
      (struct orbita_integrator *) malloc(sizeof(struct orbita_integrator));

  intgr -> delta_t = dt,
  intgr -> omega   = 0.;

  intgr -> is_adapt_dt = 0,
  intgr -> eps         = 0.;

  intgr -> order = 4;

  intgr -> wspace =
      malloc(sizeof(struct orbita_integrator_runge_kutta_4_workspace));

  intgr -> pace = & orbita_integrator_runge_kutta_4_step;

  intgr -> kill = & orbita_integrator_runge_kutta_4_kill;

  return intgr;
}

struct orbita_integrator *
orbita_integrator_runge_kutta_4_rot_init(double dt, double omega)
{
  struct orbita_integrator * intgr =
      (struct orbita_integrator *) malloc(sizeof(struct orbita_integrator));

  intgr -> delta_t = dt,
  intgr -> omega   = omega;

  intgr -> is_adapt_dt = 0,
  intgr -> eps         = 0.;

  intgr -> order = 4;

  intgr -> wspace =
      malloc(sizeof(struct orbita_integrator_runge_kutta_4_workspace));

  intgr -> pace = & orbita_integrator_runge_kutta_4_rot_step;

  intgr -> kill = & orbita_integrator_runge_kutta_4_kill;

  return intgr;
}


