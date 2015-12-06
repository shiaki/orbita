
#include "stdlib.h"
#include "stdio.h"
#include "stddef.h"
#include "math.h"

#include "utils.h"
#include "potential.h"
#include "sphintp.h"

#include "assert.h"

/*
  Locating the point on an axis
*/
int
_sphintp_axbsch(int N, double * X, double v)
{
  // indices for the left/right/middle point
  int i = 0, j = N - 1;
  int k = (i + j) / 2;

  // out-of-range check
  if(v <  X[i]) return 0;
  if(v >= X[j]) return N;

  // start binary search
  while(j - i > 1)
  {
    if(v <  X[k]) j = k;
    else          i = k;
    k = (i + j) / 2;
  }

  // return the right point idx
  return j;
}

double
_sphintp_trilinear(double * data, double fill,
           int      nI,   int      nJ,  int      nK,
           double * axI,  double * axJ, double * axK, // mono increasing
           double   pI,   double   pJ,  double   pK)
{
  int
    iI = _sphintp_axbsch(nI, axI, pI),
    iJ = _sphintp_axbsch(nJ, axJ, pJ),
    iK = _sphintp_axbsch(nK, axK, pK);

  if(iI == 0  | iJ == 0  | iK == 0 |
     iI == nI | iJ == nJ | iK == nK) return fill;

  double
    aI = axI[iI - 1], aJ = axJ[iJ - 1], aK = axK[iK - 1],
    bI = axI[iI],     bJ = axJ[iJ],     bK = axK[iK]    ;

  double
    p0 = data[(iI - 1) * nJ * nK + (iJ - 1) * nK + (iK - 1)],
    p1 = data[(iI - 1) * nJ * nK + (iJ - 1) * nK + (iK    )],
    p2 = data[(iI - 1) * nJ * nK + (iJ    ) * nK + (iK - 1)],
    p3 = data[(iI - 1) * nJ * nK + (iJ    ) * nK + (iK    )],
    p4 = data[(iI    ) * nJ * nK + (iJ - 1) * nK + (iK - 1)],
    p5 = data[(iI    ) * nJ * nK + (iJ - 1) * nK + (iK    )],
    p6 = data[(iI    ) * nJ * nK + (iJ    ) * nK + (iK - 1)],
    p7 = data[(iI    ) * nJ * nK + (iJ    ) * nK + (iK    )];

  double
    dI = (pI - aI) / (bI - aI),
    dJ = (pJ - aJ) / (bJ - aJ),
    dK = (pK - aK) / (bK - aK);

  double
    v0 = (1. - dI) * (1. - dJ) * (1. - dK),
    v1 = (1. - dI) * (1. - dJ) * (     dK),
    v2 = (1. - dI) * (     dJ) * (1. - dK),
    v3 = (1. - dI) * (     dJ) * (     dK),
    v4 = (     dI) * (1. - dJ) * (1. - dK),
    v5 = (     dI) * (1. - dJ) * (     dK),
    v6 = (     dI) * (     dJ) * (1. - dK),
    v7 = (     dI) * (     dJ) * (     dK);

  return
    (v0 * p0 + v1 * p1 + v2 * p2 + v3 * p3 +
     v4 * p4 + v5 * p5 + v6 * p6 + v7 * p7);
}



struct _orbita_models_sphintp_trilinear_paramset *
_orbita_models_sphintp_trilinear_readdump(const char * filename,
                                          const int    N_r_pts,
                                          const int    N_theta_pts,
                                          const int    N_phi_pts,
                                          const double r_max)
{
  int ir, itheta, iphi;

  FILE * dump_file = fopen(filename, "rb");
  if(!dump_file) {orbita_err("potential file does not exist."); return NULL;}

  // make workspace
  struct _orbita_models_sphintp_trilinear_paramset * ws =
    (struct _orbita_models_sphintp_trilinear_paramset *)
    malloc(sizeof(struct _orbita_models_sphintp_trilinear_paramset));

  ws -> r_pts     = N_r_pts,
  ws -> theta_pts = N_theta_pts,
  ws -> phi_pts   = N_phi_pts;

  ws -> r_ax     = (double *) malloc(sizeof(double) * N_r_pts),
  ws -> theta_ax = (double *) malloc(sizeof(double) * N_theta_pts),
  ws -> phi_ax   = (double *) malloc(sizeof(double) * N_phi_pts);

  //ws -> dphi   = fmod(delta_phi, Pi);

  // radial axis points
  double alpha_R = log(r_max + 1.) / (N_r_pts - 1);
  for(ir = 0; ir < N_r_pts; ++ ir)
    (ws -> r_ax)[ir] = exp(ir * alpha_R) - 1.;

  // theta points
  for(itheta = 0; itheta < N_theta_pts; ++ itheta)
    (ws -> theta_ax)[itheta] = itheta * Pi / (N_theta_pts - 1);

  // azimuthal axis points
  for(iphi = 0; iphi < N_phi_pts; ++ iphi)
    (ws -> phi_ax)[iphi] = iphi * (2. * Pi) / (N_phi_pts - 1);//,

  // allocate array space
  size_t ws_gridsize = N_r_pts * N_theta_pts * N_phi_pts;

  ws -> Fx_val  = (double *) malloc(sizeof(double) * ws_gridsize),
  ws -> Fy_val  = (double *) malloc(sizeof(double) * ws_gridsize),
  ws -> Fz_val  = (double *) malloc(sizeof(double) * ws_gridsize),
  ws -> psi_val = (double *) malloc(sizeof(double) * ws_gridsize);

  // prepare temp array for i/o
  //double * arr_t = (double *) malloc(sizeof(double) * ws_gridsize);

  // Now read file

  // read psi into memory and copy to workspace
  size_t rec_size;
  rec_size = fread((void *)(ws -> psi_val), sizeof(double), ws_gridsize, dump_file);
  rec_size = fread((void *)(ws -> Fx_val ), sizeof(double), ws_gridsize, dump_file);
  rec_size = fread((void *)(ws -> Fy_val ), sizeof(double), ws_gridsize, dump_file);
  rec_size = fread((void *)(ws -> Fz_val ), sizeof(double), ws_gridsize, dump_file);

  assert(rec_size == ws_gridsize);

  fclose(dump_file);

  return ws;
}

double
_orbita_models_sphintp_trilinear_rho(const void * par, const double * pos)
{
  printf("sphintp_trilinear does not support rho!\n");
  return -1e6;
}

double
_orbita_models_sphintp_trilinear_Mr(const void * par, const double R)
{
  printf("sphintp_trilinear does not support Mr!\n");
  return -1e6;
}

double
_orbita_models_sphintp_trilinear_Mt(const void * par)
{
  printf("sphintp_trilinear does not support Mt!\n");
  return -1e6;
}

double
_orbita_models_sphintp_trilinear_psi(const void * par, const double * pos)
{

  struct _orbita_models_sphintp_trilinear_paramset * w =
    (struct _orbita_models_sphintp_trilinear_paramset *) par;

  //double dphi = w -> dphi;

  double r     = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]),
         theta = acos(pos[2] / r),
         phi   = fmod(atan2(pos[1], pos[0]) + 2. * Pi, 2. * Pi);

  return -_sphintp_trilinear(w -> psi_val, 0.,
                     w -> r_pts, w -> theta_pts, w -> phi_pts,
                     w -> r_ax,  w -> theta_ax,  w -> phi_ax,
                          r,          theta,          phi);
}

int
_orbita_models_sphintp_trilinear_f(const void * par,
                                  const double * pos, double * F)
{
  struct _orbita_models_sphintp_trilinear_paramset * w =
    (struct _orbita_models_sphintp_trilinear_paramset *) par;

  double r     = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]),
         theta = acos(pos[2] / r),
         phi   = fmod(atan2(pos[1], pos[0]) + 2. * Pi, 2. * Pi);

  double Fx   = _sphintp_trilinear(w -> Fx_val, 0.,
                           w -> r_pts, w -> theta_pts, w -> phi_pts,
                           w -> r_ax,  w -> theta_ax,  w -> phi_ax,
                                r,          theta,          phi),
         Fy   = _sphintp_trilinear(w -> Fy_val, 0.,
                           w -> r_pts, w -> theta_pts, w -> phi_pts,
                           w -> r_ax,  w -> theta_ax,  w -> phi_ax,
                                r,          theta,          phi),
         Fz = _sphintp_trilinear(w -> Fz_val, 0.,
                           w -> r_pts, w -> theta_pts, w -> phi_pts,
                           w -> r_ax,  w -> theta_ax,  w -> phi_ax,
                                r,          theta,          phi);

  F[0] = Fx, F[1] = Fy, F[2] = Fz;

  return 0;
}

int
_orbita_models_sphintp_kill(struct orbita_potential * psi)
{
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> r_ax),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> theta_ax),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> phi_ax),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> Fx_val),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> Fy_val),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> Fz_val),
  free(((struct _orbita_models_sphintp_trilinear_paramset *)(psi -> param)) -> psi_val),
  free(psi -> param);

  return 0;
}

struct orbita_potential *
orbita_potential_sphintp_trilinear(const char *  dump_filename,
                                   const int     r_pts,
                                   const int     theta_pts,
                                   const int     phi_pts,
                                   const double  r_max)
{
  struct orbita_potential * psi = (struct orbita_potential *)
      malloc(sizeof(struct orbita_potential));

  psi -> param = (void *) _orbita_models_sphintp_trilinear_readdump
      (dump_filename, r_pts, theta_pts, phi_pts, r_max);

  psi -> symmetry     = 0,
  psi -> coordinates  = ORBITA_COORD_SPHERICAL;

  psi -> rho = & _orbita_models_sphintp_trilinear_rho,
  psi -> psi = & _orbita_models_sphintp_trilinear_psi,
  psi -> f   = & _orbita_models_sphintp_trilinear_f,
  psi -> Mr  = & _orbita_models_sphintp_trilinear_Mr,
  psi -> Mt  = & _orbita_models_sphintp_trilinear_Mt;

  psi -> is_variable = 0,
  psi -> init_param  = NULL,
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_models_sphintp_kill;

  psi -> is_composite = 0;

  return psi;
}
