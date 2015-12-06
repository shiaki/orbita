
#include "stdlib.h"
#include "stdio.h"
#include "stddef.h"
#include "math.h"

#include "utils.h"
#include "potential.h"
#include "fcsintp.h"

#include "assert.h"

/* within the header */

int
_axbsch(int N, float * X, double v)
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

int
_make_ghost_layer(float * src, float * grid,
                  int N_R_pts, int N_Z_pts, int N_phi_pts)
{
  int ir, iz, iphi;

  for(ir = 0; ir < N_R_pts; ++ ir)
    for(iz = 0; iz < N_Z_pts; ++ iz)
      for(iphi = 0; iphi <= N_phi_pts; ++ iphi)
        grid[ ir * (N_Z_pts * (N_phi_pts + 1))
             + iz * (N_phi_pts + 1) + iphi]
        = src[ir * (N_Z_pts * N_phi_pts)
             + iz *  N_phi_pts +     (iphi % N_phi_pts)];

  return 0;
}

struct _orbita_models_fcsp3d_trilinear_paramset *
_orbita_models_fcsp3d_trilinear_readdump(const char * dump_filename,
                                   const int    N_R_pts,
                                   const int    N_Z_pts,
                                   const int    N_phi_pts,
                                   const float  R_scale,
                                   const float  Z_scale,
                                   const float  delta_phi)
{
  int ir, iz, iphi;

  FILE * dump_file = fopen(dump_filename, "rb");
  if(!dump_file) return NULL;

  // make workspace
  struct _orbita_models_fcsp3d_trilinear_paramset * ws =
      (struct _orbita_models_fcsp3d_trilinear_paramset *)
      malloc(sizeof(struct _orbita_models_fcsp3d_trilinear_paramset));

  ws -> R_pts   = N_R_pts,
  ws -> Z_pts   = N_Z_pts,
  ws -> phi_pts = N_phi_pts + 1;

  ws -> R_ax    = (float *) malloc(sizeof(float) * N_R_pts),
  ws -> Z_ax    = (float *) malloc(sizeof(float) * N_Z_pts),
  ws -> phi_ax  = (float *) malloc(sizeof(float) * (N_phi_pts + 1));

  ws -> dphi    = fmod(delta_phi, Pi);

  // radial axis points
  float alpha_R = 2. * Pi / N_phi_pts;
  for(ir = 0; ir < N_R_pts; ++ ir)
    (ws -> R_ax)[ir] = (exp(ir * alpha_R) - 1.) / R_scale;

  // vertical axis points
  for(iz = 0; iz < N_Z_pts; ++ iz)
    (ws -> Z_ax)[iz] = (iz - (N_Z_pts - 1.) / 2.) / Z_scale;

  // azimuthal axis points
  for(iphi = 0; iphi <= N_phi_pts; ++ iphi)
    (ws -> phi_ax)[iphi] = iphi * (2. * Pi) / N_phi_pts;//,
    //printf("%u %f\n", iphi, (ws -> phi_ax)[iphi]);

  // allocate array space
  size_t fcs_gridsize = N_R_pts * N_Z_pts * N_phi_pts,
         ws_gridsize  = N_R_pts * N_Z_pts * (N_phi_pts + 1);

  ws -> FR_val   = (float *) malloc(sizeof(float) * ws_gridsize),
  ws -> FZ_val   = (float *) malloc(sizeof(float) * ws_gridsize),
  ws -> Fphi_val = (float *) malloc(sizeof(float) * ws_gridsize),
  ws -> psi_val  = (float *) malloc(sizeof(float) * ws_gridsize);

  // prepare temp array for i/o
  float * arr_t = (float *) malloc(sizeof(float) * fcs_gridsize);

  // Now read fcs file.

  // used to skip fortran record mark
  unsigned int rec_mark; size_t rec_size;

  // read FR into memory and copy to workspace
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);
  rec_size = fread(arr_t, sizeof(float), fcs_gridsize, dump_file);
  _make_ghost_layer(arr_t, ws -> FR_val, N_R_pts, N_Z_pts, N_phi_pts);
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);

  // read Fphi into memory and copy to workspace
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);
  rec_size = fread(arr_t, sizeof(float), fcs_gridsize, dump_file);
  _make_ghost_layer(arr_t, ws -> Fphi_val, N_R_pts, N_Z_pts, N_phi_pts);
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);

  // read FZ into memory and copy to workspace
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);
  rec_size = fread(arr_t, sizeof(float), fcs_gridsize, dump_file);
  _make_ghost_layer(arr_t, ws -> FZ_val, N_R_pts, N_Z_pts, N_phi_pts);
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);

  // read psi into memory and copy to workspace
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);
  rec_size = fread(arr_t, sizeof(float), fcs_gridsize, dump_file);
  _make_ghost_layer(arr_t, ws -> psi_val, N_R_pts, N_Z_pts, N_phi_pts);
  fread(& rec_mark, sizeof(unsigned int), 1, dump_file);

  fclose(dump_file);

  return ws;
}

float
_trilinear(float * data, float fill,
           int      nI, int      nJ, int      nK,
           float * axI, float * axJ, float * axK, // mono increasing
           double   pI, double   pJ, double   pK)
{
  int
    iI = _axbsch(nI, axI, pI),
    iJ = _axbsch(nJ, axJ, pJ),
    iK = _axbsch(nK, axK, pK);

  if((iI == 0 ) | (iJ == 0 ) | (iK == 0 )|
     (iI == nI) | (iJ == nJ) | (iK == nK)) return fill;

  float
    aI = axI[iI - 1], aJ = axJ[iJ - 1], aK = axK[iK - 1],
    bI = axI[iI],     bJ = axJ[iJ],     bK = axK[iK]    ;

  float
    p0 = data[(iI - 1) * nJ * nK + (iJ - 1) * nK + (iK - 1)],
    p1 = data[(iI - 1) * nJ * nK + (iJ - 1) * nK + (iK    )],
    p2 = data[(iI - 1) * nJ * nK + (iJ    ) * nK + (iK - 1)],
    p3 = data[(iI - 1) * nJ * nK + (iJ    ) * nK + (iK    )],
    p4 = data[(iI    ) * nJ * nK + (iJ - 1) * nK + (iK - 1)],
    p5 = data[(iI    ) * nJ * nK + (iJ - 1) * nK + (iK    )],
    p6 = data[(iI    ) * nJ * nK + (iJ    ) * nK + (iK - 1)],
    p7 = data[(iI    ) * nJ * nK + (iJ    ) * nK + (iK    )];

  float
    dI = (pI - aI) / (bI - aI),
    dJ = (pJ - aJ) / (bJ - aJ),
    dK = (pK - aK) / (bK - aK);

  float
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

double
_orbita_models_fcsp3d_trilinear_rho(const void * par, const double * pos)
{
  printf("fcsp3d_trilinear does not support rho!\n");
  return -1e6;
}

double
_orbita_models_fcsp3d_trilinear_Mr(const void * par, const double R)
{
  printf("fcsp3d_trilinear does not support Mr!\n");
  return -1e6;
}

double
_orbita_models_fcsp3d_trilinear_Mt(const void * par)
{
  printf("fcsp3d_trilinear does not support Mt!\n");
  return -1e6;
}

double
_orbita_models_fcsp3d_trilinear_psi(const void * par, const double * pos)
{
  struct _orbita_models_fcsp3d_trilinear_paramset * w =
    (struct _orbita_models_fcsp3d_trilinear_paramset *) par;

  double dphi = w -> dphi;

  double R    = sqrt(pos[0] * pos[0] + pos[1] * pos[1]),
         phi  = fmod(atan2(pos[1], pos[0]) - dphi + 2.* Pi, 2. * Pi),
         Z    = pos[2];

  return -_trilinear(w -> psi_val, 0.,
                     w -> R_pts, w -> Z_pts, w -> phi_pts,
                     w -> R_ax,  w -> Z_ax,  w -> phi_ax,
                          R,          Z,          phi);
}

int
_orbita_models_fcsp3d_trilinear_f(const void * par,
                                  const double * pos, double * F)
{
  struct _orbita_models_fcsp3d_trilinear_paramset * w =
    (struct _orbita_models_fcsp3d_trilinear_paramset *) par;

  double dphi = w -> dphi;

  double R    = sqrt(pos[0] * pos[0] + pos[1] * pos[1]),
         phi  = fmod(atan2(pos[1], pos[0]) - dphi + 2.* Pi, 2. * Pi),
         Z    = pos[2];

  double FR   = _trilinear(w -> FR_val, 0.,
                           w -> R_pts, w -> Z_pts, w -> phi_pts,
                           w -> R_ax,  w -> Z_ax,  w -> phi_ax,
                                R,          Z,          phi),
         FZ   = _trilinear(w -> FZ_val, 0.,
                           w -> R_pts, w -> Z_pts, w -> phi_pts,
                           w -> R_ax,  w -> Z_ax,  w -> phi_ax,
                                R,          Z,          phi),
         Fphi = _trilinear(w -> Fphi_val, 0.,
                           w -> R_pts, w -> Z_pts, w -> phi_pts,
                           w -> R_ax,  w -> Z_ax,  w -> phi_ax,
                                R,          Z,          phi);

  double cos_phi = pos[0] / R,
         sin_phi = pos[1] / R;

  F[0] = FR * cos_phi - Fphi * sin_phi,
  F[1] = FR * sin_phi + Fphi * cos_phi,
  F[2] = FZ;

  return 0;
}

int
_orbita_models_fcsp3d_kill(struct orbita_potential * psi)
{
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> R_ax),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> Z_ax),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> phi_ax),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> FR_val),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> FZ_val),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> Fphi_val),
  free(((struct _orbita_models_fcsp3d_trilinear_paramset *)
        (psi -> param)) -> psi_val),
  free(psi -> param);
  //free(psi);
  return 0;
}

struct orbita_potential *
orbita_potential_fcsp3d_trilinear(const char *  fcs_filename,
                                  const int     R_pts,
                                  const int     Z_pts,
                                  const int     phi_pts,
                                  const double  R_scale,
                                  const double  Z_scale,
                                  const double  delta_phi)
{
  struct orbita_potential * psi = (struct orbita_potential *)
      malloc(sizeof(struct orbita_potential));

  psi -> param = (void *) _orbita_models_fcsp3d_trilinear_readdump
      (fcs_filename, R_pts, Z_pts, phi_pts, R_scale, Z_scale, delta_phi);

  psi -> symmetry     = 0,
  psi -> coordinates  = ORBITA_COORD_CYLINDRICAL;

  psi -> rho = & _orbita_models_fcsp3d_trilinear_rho,
  psi -> psi = & _orbita_models_fcsp3d_trilinear_psi,
  psi -> f   = & _orbita_models_fcsp3d_trilinear_f,
  psi -> Mr  = & _orbita_models_fcsp3d_trilinear_Mr,
  psi -> Mt  = & _orbita_models_fcsp3d_trilinear_Mt;

  psi -> is_variable = 0,
  psi -> init_param  = NULL,
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_models_fcsp3d_kill;

  psi -> is_composite = 0;

  return psi;
}


//int
//main(int argc, char ** argv)
//{
  // Test binary search
  /*
  double a[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
  int i = _axbsch(8, a, 4.);
  printf("i = %u\n", i);
  */

  // Test loading file
  /*
    _orbita_models_fcsp3d_trilinear_readdump(const char * dump_filename,
                                       const int    N_R_pts,
                                       const int    N_Z_pts,
                                       const int    N_phi_pts,
                                       const float  R_scale,
                                       const float  Z_scale)
  */
  /*
  struct _orbita_models_fcsp3d_trilinear_paramset * ws =
    _orbita_models_fcsp3d_trilinear_readdump("../run3042.fcs0800",
        58, 375, 64, 10., 100.);

  printf("%u\n", ws -> phi_pts);
  printf("%f\n", ws -> psi_val[127]);
  */

  //

  /*
  struct _orbita_models_fcsp3d_trilinear_paramset * ws =
    _orbita_models_fcsp3d_trilinear_readdump("../run3042.fcs0800",
        58, 375, 64, 10., 100.);

  double pos_t[3] = {2., 0., 0.};
  double psi_t = _orbita_models_fcsp3d_trilinear_psi(ws, pos_t);

  printf("%f\n", psi_t);

  int idx = _axbsch(65, ws -> phi_ax, 1.5);
  printf("Located: %u\n", idx);
  */

  /*
  double f = fmod(1.507, 2. * Pi);
  printf("f = %f\n", f);
  */

  /*
  struct _orbita_models_fcsp3d_trilinear_paramset * ws =
    _orbita_models_fcsp3d_trilinear_readdump("../run3042.fcs0800",
        58, 375, 64, 10., 100.);

  #define NPTS 1025
  double ax[NPTS], pos[3], val; int i, j;
  for(i = 0; i < NPTS; ++ i) ax[i] = -2. + (4. / (NPTS - 1)) * i;

  for(i = 0; i < NPTS; ++ i)
  {
    for(j = 0; j < NPTS; ++ j)
    {
      pos[0] = ax[i], pos[1] = ax[j], pos[2] = 0.1;
      //printf("pos: %f, %f, %f\n", pos[0], pos[1], pos[2]);
      val = _orbita_models_fcsp3d_trilinear_psi(ws, pos);
      printf("%f ", val);
    }
    printf("\n");
  }
  */

  /*
  for(j = 10; j < 15; ++ j)
  {
    //pos[0] = cos(j * Pi / 2. / NPTS),
    //pos[1] = sin(j * Pi / 2. / NPTS),
    pos[0] = 5. + (float)j / NPTS,
    pos[1] = 5.01 + (float)j / NPTS,
    pos[2] = 0.0213;
    //printf("pos: %f, %f, %f\n", pos[0], pos[1], pos[2]);
    val = _orbita_models_fcsp3d_trilinear_psi(ws, pos);
    printf("%f ", val);
  }
  */
  // dump a page from psi array
  //float m; int ip, iz, ir;
  /*
  for(ip = 0; ip < 65; ++ ip)
  {
    for(iz = 0; iz < 375; ++ iz)
    {
      m = (ws -> psi_val)[20 * (65 * 375) + iz * 65 + ip];
      printf("%f ", m);
    }
    printf("\n");
  }
  */
  /*
  for(ir = 0; ir < 58; ++ ir)
  {
    for(iz = 0; iz < 375; ++ iz)
    {
      m = (ws -> psi_val)[ir * (65 * 375) + iz * 65 + 12];
      printf("%f ", m);
    }
    printf("\n");
  }
  */
  /*
  for(ip = 0; ip < 65; ++ ip)
  {
    for(ir = 0; ir < 58; ++ ir)
    {
      m = (ws -> psi_val)[ir * (65 * 375) + 187 * 65 + ip];
      printf("%f ", m);
    }
    printf("\n");
  }
  */
//}
