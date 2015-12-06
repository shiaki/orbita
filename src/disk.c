
/*
  disk.c
    Models of for galaxy disks.
*/

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "potential.h"

#include "disk.h"

#define Sq(X) ((X) * (X))

/*=========================================================
  I. Miyamoto-Nagai Disk;
    psi = -G M / sqrt(x^2 + y^2 + (A + sqrt(B^2 + z^2))^2)
    Fx  = -(G M x)/(A^2+2 A sqrt(B^2+z^2)+B^2+x^2+y^2+z^2)^(3/2)
    Fy  = -(G M y)/(A^2+2 A sqrt(B^2+z^2)+B^2+x^2+y^2+z^2)^(3/2)
    Fz  = -(G M z)/V - (A G M z)/(sqrt(B^2+z^2) V)
=========================================================*/

int
_orbita_potential_miyamoto_disk_f(const void * par, const double * pos, double * F)
{
  double M = ((double *) par)[0], A = ((double *) par)[1], B_sq = ((double *) par)[2],
         x_sq = pos[0] * pos[0],  y_sq = pos[1] * pos[1],  z_sq = pos[2] * pos[2];

  double W = sqrt(B_sq + z_sq),
         U = pow(Sq(A) + 2. * A * W + B_sq + x_sq + y_sq + z_sq, 1.5),
         V = pow(Sq(A + W) + x_sq + y_sq, 1.5);

  F[0] = -(M * pos[0]) / U,
  F[1] = -(M * pos[1]) / U,
  F[2] = -(M * pos[2]) / V - (A * M * pos[2]) / (W * V);

  return 0;
}

double
_orbita_potential_miyamoto_disk_psi(const void * par, const double * pos)
{
  double M = ((double *) par)[0], A = ((double *) par)[1], B_sq = ((double *) par)[2],
         x_sq = pos[0] * pos[0],  y_sq = pos[1] * pos[1],  z_sq = pos[2] * pos[2];
  return -M / sqrt(x_sq + y_sq + Sq(A + sqrt(B_sq + z_sq)));
}

double
_orbita_potential_miyamoto_disk_rho(const void * par, const double * pos)
{
  orbita_err("Not implemented");
  return 0. / 0.;
}

double
_orbita_potential_miyamoto_disk_Mr(const void * par, const double R)
{
  orbita_err("M(R) not implemented for Miyamoto-Nagai disk.");
  return -1;
}

double
_orbita_potential_miyamoto_disk_Mt(const void * par)
{
  return ((double *) par)[0];
}

int
_orbita_potential_miyamoto_disk_kill(struct orbita_potential * psi)
{
  free(psi -> param); return 0;
}

struct orbita_potential *
orbita_potential_miyamoto_disk(double M, double A, double B)
{
  struct orbita_potential * psi =
    (struct orbita_potential *) malloc(sizeof(struct orbita_potential));
  double * ws = (double *) malloc(sizeof(double) * 3);

  ws[0] = M, ws[1] = A, ws[2] = B * B;
  psi -> param = (void *) ws;

  // set other parameters
  psi -> symmetry     = 0,
  psi -> coordinates  = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_potential_miyamoto_disk_rho,
  psi -> psi = & _orbita_potential_miyamoto_disk_psi,
  psi -> f   = & _orbita_potential_miyamoto_disk_f,
  psi -> Mr  = & _orbita_potential_miyamoto_disk_Mr,
  psi -> Mt  = & _orbita_potential_miyamoto_disk_Mt;

  psi -> is_variable = 0,
  psi -> init_param  = NULL,
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_potential_miyamoto_disk_kill;

  psi -> is_composite = 0;

  return psi;
}
