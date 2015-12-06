
/*
  ferrers.c
    Ferrers bar model.
*/

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "potential.h"

#include <gsl/gsl_sf.h>

#include "ferrers.h"

#define Sq(X) ((X) * (X))

double
_lambda(double x_sq, double y_sq, double z_sq, double a_sq, double b_sq, double c_sq)
{
  #define EQS(u) (x_sq / (a_sq + u) + y_sq / (b_sq + u) + z_sq / (c_sq + u) - 1.)
  double J = 0, K = 0., L = (x_sq + y_sq + z_sq);
  for(;EQS(L) > 0;) L *= 2.;//, printf("doubled!\n");

  // solve the eqs using bisection method.
  for(;L - J > 1.e-8;)
    {
      K = (L + J) / 2.;
      //printf("%f, %f, %f\n", J, K, L);
      if(EQS(K) > 0) J = K;
      else L = K;
    }

  return K;

  #undef EQS
}

int
_orbita_potential_ferrers2_evaluate_wijk(double x, double y, double z,
                                         double a, double b, double c,
                                         double * W, int init,
                                         struct orbita_potential_ferrers2_paramset * ws)
{
  double a_sq = a * a, b_sq = b * b, c_sq = c * c,
         x_sq = x * x, y_sq = y * y, z_sq = z * z;

  // outside the bar region?
  if(x_sq / a_sq + y_sq / b_sq + z_sq / c_sq >= 1 || init)
    {
      // using the root-finding eqs. this is NOT the unique positive root.
      // now use bisection method instead.
      /*
      double A = (a_sq + b_sq + c_sq) - (x_sq + y_sq + z_sq),
             B = (a_sq * b_sq + b_sq * c_sq + c_sq * a_sq) -
                 (x_sq * (b_sq + c_sq) + y_sq * (c_sq + a_sq) + z_sq * (a_sq + b_sq)),
             C = (a_sq * b_sq * c_sq) -
                 (x_sq * b_sq * c_sq + a_sq * y_sq * c_sq + a_sq * b_sq * z_sq);

      double A_sq = A * A, B_sq = B * B, C_sq = C * C;
      double A_cb = A_sq * A, B_cb = B_sq * B;
      double CBR_2 = pow(2., 1. / 3.);

      double lambda =
        pow(-2. * A_cb + 3. * sqrt(3.) * sqrt(4. * A_cb * C - A_sq * B_sq - 18. * A * B * C
          + 4. * B_cb + 27. * C_sq) + 9. * A * B - 27. * C, 1. / 3.) / (3. * CBR_2)
          - (CBR_2 * (2187. * B - 729. * A_sq)) / (2187. * pow(-2. * A_cb + 3. * sqrt(3.)
          * sqrt(4. * A_cb * C - A_sq * B_sq - 18. * A * B * C + 4. * B_cb + 27. * C_sq)
          + 9. * A * B - 27. * C, 1. / 3.)) - A / 3.;
      */

      double lambda = _lambda(x_sq, y_sq, z_sq, a_sq, b_sq, c_sq);
      if(init) lambda = 0.;

      double p = asin(sqrt((a_sq - c_sq) / (a_sq + lambda))),
             k = sqrt((a_sq - b_sq) / (a_sq - c_sq));

      /*
      double F = gsl_sf_ellint_F(p, k, GSL_PREC_DOUBLE),
             E = gsl_sf_ellint_E(p, k, GSL_PREC_DOUBLE);
      */

      double F = gsl_sf_ellint_F(p, k, GSL_PREC_SINGLE),
             E = gsl_sf_ellint_E(p, k, GSL_PREC_SINGLE);

      double delta_lambda = sqrt((a_sq + lambda) * (b_sq + lambda) * (c_sq + lambda));

      /* W_000 */ W[ 0] = 2. * F / sqrt(a_sq - c_sq),
      /* W_100 */ W[ 1] = 2. * (F - E) / ((a_sq - b_sq) * sqrt(a_sq - c_sq)),
      /* W_001 */ W[ 3] = (2. / (b_sq - c_sq)) * sqrt((b_sq + lambda) / ((a_sq + lambda) * (c_sq + lambda)))
                          -2. * E / ((b_sq - c_sq) * sqrt(a_sq - c_sq)),
      /* W_010 */ W[ 2] = 2. / delta_lambda - W[1] - W[3],
      /* W_110 */ W[ 4] = (W[2] - W[1]) / (a_sq - b_sq),
      /* W_011 */ W[ 5] = (W[3] - W[2]) / (b_sq - c_sq),
      /* W_101 */ W[ 6] = (W[1] - W[3]) / (c_sq - a_sq),
      /* W_200 */ W[ 7] = (2. / (delta_lambda * (a_sq + lambda)) - W[4] - W[6]) / 3.,
      /* W_020 */ W[ 8] = (2. / (delta_lambda * (b_sq + lambda)) - W[5] - W[4]) / 3.,
      /* W_002 */ W[ 9] = (2. / (delta_lambda * (c_sq + lambda)) - W[6] - W[5]) / 3.,
      /* W_111 */ W[10] = (W[4] - W[5]) / (c_sq - a_sq),
      /* W_120 */ W[11] = (W[8] - W[4]) / (a_sq - b_sq),
      /* W_012 */ W[12] = (W[9] - W[5]) / (b_sq - c_sq),
      /* W_201 */ W[13] = (W[7] - W[6]) / (c_sq - a_sq),
      /* W_210 */ W[14] = (W[7] - W[4]) / (b_sq - a_sq), // NOTE Fixed. Nov. 15, 2015
      /* W_021 */ W[15] = (W[8] - W[5]) / (c_sq - b_sq), //
      /* W_102 */ W[16] = (W[9] - W[6]) / (a_sq - c_sq), //
      /* W_300 */ W[17] = (2. / (delta_lambda * Sq(a_sq + lambda)) - W[14] - W[13]) / 5.,
      /* W_030 */ W[18] = (2. / (delta_lambda * Sq(b_sq + lambda)) - W[15] - W[11]) / 5.,
      /* W_003 */ W[19] = (2. / (delta_lambda * Sq(c_sq + lambda)) - W[16] - W[12]) / 5.;
    }
  else
    {
      int i_ijk;
      for(i_ijk = 0; i_ijk < 20; ++ i_ijk)
        W[i_ijk] = (ws -> W_i)[i_ijk];
    }

  return 0;
}

double
_orbita_potential_ferrers2_psi(const void * par, const double * pos)
{
  struct orbita_potential_ferrers2_paramset * ws =
    (struct orbita_potential_ferrers2_paramset *) par;

  double W[20];
  _orbita_potential_ferrers2_evaluate_wijk(
      pos[0], pos[1], pos[2], ws -> a, ws -> b, ws -> c, W, 0, ws);

  double x_sq = pos[0] * pos[0], y_sq = pos[1] * pos[1], z_sq = pos[2] * pos[2];

  /*
  return (-(ws -> C) / 6.) *
         (
             W[0] - 6. * x_sq * y_sq * z_sq * W[10] 
           + x_sq * (x_sq * (3. * W[7] - x_sq * W[17]) 
           + 3. * (y_sq * (2 * W[4] - y_sq * W[11] - x_sq * W[14]) - W[1])) 
           + y_sq * (y_sq * (3. * W[8] - y_sq * W[18])
           + 3. * (z_sq * (2 * W[5] - z_sq * W[12] - y_sq * W[15]) - W[2])) 
           + z_sq * (z_sq * (3. * W[9] - z_sq * W[19]) 
           + 3. * (x_sq * (2 * W[6] - x_sq * W[13] - z_sq * W[16]) - W[3])) 
         );
  */
  
  // This is to test if Del^2 Psi == 4 Pi G rho. Returns rho from W_jkl.
  return ((ws -> a) * (ws -> b) * (ws -> c) * (ws -> rho) / 2.) *
          (
                  (W[1] + W[2] + W[3])
            - 2. * ((y_sq + x_sq) * W[4] + (z_sq + y_sq) * W[5] + (x_sq + z_sq) * W[6])
            - 6. * (x_sq * W[7] + y_sq * W[8] + z_sq * W[9])
            + 2. * (x_sq * y_sq + y_sq * z_sq + z_sq * x_sq) * W[10]
            + (y_sq * y_sq + 6. * x_sq * y_sq) * W[11]
            + (z_sq * z_sq + 6. * y_sq * z_sq) * W[12]
            + (x_sq * x_sq + 6. * z_sq * x_sq) * W[13]
            + (6. * x_sq * y_sq + x_sq * x_sq) * W[14]
            + (6. * y_sq * z_sq + y_sq * y_sq) * W[15]
            + (6. * z_sq * x_sq + z_sq * z_sq) * W[16]
            + 5. * (x_sq * x_sq * W[17] + y_sq * y_sq * W[18] + z_sq * z_sq * W[19])
          ); 
}

int
_orbita_potential_ferrers2_f(const void * par,
                             const double * pos, double * F)
{
  struct orbita_potential_ferrers2_paramset * ws =
    (struct orbita_potential_ferrers2_paramset *) par;

  double W[20];
  _orbita_potential_ferrers2_evaluate_wijk(
    pos[0], pos[1], pos[2], ws -> a, ws -> b, ws -> c, W, 0, ws);

  double x_sq = pos[0] * pos[0], y_sq = pos[1] * pos[1], z_sq = pos[2] * pos[2]; //

  F[0] = -pos[0] * (ws -> C) * //
          (
            W[1] + x_sq * (x_sq * W[17] + 2. * (y_sq * W[14] - W[7]))
                 + y_sq * (y_sq * W[11] + 2. * (z_sq * W[10] - W[4]))
                 + z_sq * (z_sq * W[16] + 2. * (x_sq * W[13] - W[6]))
          ),
  F[1] = -pos[1] * (ws -> C) * //
          (
            W[2] + x_sq * (x_sq * W[14] + 2. * (y_sq * W[11] - W[4]))
                 + y_sq * (y_sq * W[18] + 2. * (z_sq * W[15] - W[8]))
                 + z_sq * (z_sq * W[12] + 2. * (x_sq * W[10] - W[5]))
          ),
  F[2] = -pos[2] * (ws -> C) * //
          (
            W[3] + x_sq * (x_sq * W[13] + 2. * (y_sq * W[10] - W[6]))
                 + y_sq * (y_sq * W[15] + 2. * (z_sq * W[12] - W[5]))
                 + z_sq * (z_sq * W[19] + 2. * (x_sq * W[16] - W[9]))
          );

  return 0;
}

double
_orbita_potential_ferrers2_rho(const void * par, const double * pos)
{
  struct orbita_potential_ferrers2_paramset * ws =
    (struct orbita_potential_ferrers2_paramset *) par;

  double m_sq = Sq(pos[0]) / Sq(ws -> a) + Sq(pos[1]) / Sq(ws -> b) + Sq(pos[2]) / Sq(ws -> c);
  if(m_sq < 1.) return (ws -> rho) * Sq(1 - m_sq);
  else return 0.;
}

double
_orbita_potential_ferrers2_Mr(const void * par, const double R)
{
  orbita_err("Method not supported.");
  return 0.;
}

double
_orbita_potential_ferrers2_Mt(const void * par)
{
  struct orbita_potential_ferrers2_paramset * ws =
    (struct orbita_potential_ferrers2_paramset *) par;
  return (ws -> rho) * (ws -> a) * (ws -> b) * (ws -> c) * (32. * Pi / 105.);
}

int
_orbita_potential_ferrers2_kill(struct orbita_potential * psi)
{
  free(psi -> param);
  return 0;
}

struct orbita_potential *
orbita_potential_ferrers2(double rho, double a, double b, double c)
{
  struct orbita_potential * psi =
    (struct orbita_potential *) malloc(sizeof(struct orbita_potential));
  struct orbita_potential_ferrers2_paramset * ws =
    (struct orbita_potential_ferrers2_paramset *)
    malloc(sizeof(struct orbita_potential_ferrers2_paramset));
  psi -> param = (void *)ws;

  // initialize w_i
  _orbita_potential_ferrers2_evaluate_wijk(0., 0., 0., a, b, c, ws -> W_i, 1, ws);
  ws -> C = 2. * Pi * rho * a * b * c;
  ws -> a = a, ws -> b = b, ws -> c = c, ws -> rho = rho;

  // set other parameters
  psi -> symmetry     = 0,
  psi -> coordinates  = ORBITA_COORD_CARTESIAN;

  psi -> rho = & _orbita_potential_ferrers2_rho,
  psi -> psi = & _orbita_potential_ferrers2_psi,
  psi -> f   = & _orbita_potential_ferrers2_f,
  psi -> Mr  = & _orbita_potential_ferrers2_Mr,
  psi -> Mt  = & _orbita_potential_ferrers2_Mt;

  psi -> is_variable = 0,
  psi -> init_param  = NULL,
  psi -> param_func  = NULL;

  psi -> kill = & _orbita_potential_ferrers2_kill;

  psi -> is_composite = 0;

  return psi;
}
