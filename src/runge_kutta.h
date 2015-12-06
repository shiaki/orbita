
#ifndef RUNGE_KUTTA_H

struct orbita_integrator_runge_kutta_4_workspace
{
  double K1[6], K2[6], K3[6], K4[6];
  double P2[6], P3[6], P4[6];
};

int orbita_integrator_runge_kutta_4_lite_step(
    struct orbita_potential * psi, double * wi, double dt, double t,
    double omega, double * wf);

int orbita_integrator_runge_kutta_4_step(void *,
    struct orbita_potential * psi, double * wi, double dt, double t,
    double omega, double * wf);

int orbita_integrator_runge_kutta_4_rot_step(void *,
    struct orbita_potential * psi, double * wi, double dt, double t,
    double omega, double * wf);

struct orbita_integrator * orbita_integrator_runge_kutta_4_init(double dt);
struct orbita_integrator * orbita_integrator_runge_kutta_4_rot_init(double dt, double omega);
int orbita_integrator_runge_kutta_4_kill(struct orbita_integrator * intgr);


#define RUNGE_KUTTA_H
#endif
