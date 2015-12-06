
#include "stdlib.h"
#include "stdio.h"
#include "analyzer.h"
#include "poincare.h"
#include "utils.h"

int
_orbita_analyzer_poincare_section_x_peek(void * wspac, double * pt,
                                         int * i_pt, const int N_pts)
{
  // unpack the workspace
  struct _orbita_analyzer_poincare_section_wspac * ws =
    (struct _orbita_analyzer_poincare_section_wspac *) wspac;

  // having enough points on the SoS? escape.
  if(ws -> I_point + 1 == ws -> N_points) return 1;

  // skip if this is the first step of an cycle
  if(*i_pt == 0) return 0;

  // end of this orbit but still needs points
  if(*i_pt + 2 == N_pts)
    {
      if(ws -> I_cycle + 1 == ws -> N_cycles) return 1;

      int i_k;
      for(i_k = 0; i_k < 6; ++ i_k)
        pt[i_k + 6] = pt[6 * (N_pts - 2) + i_k];
      * i_pt = 0, ++ (ws -> I_cycle);

      return 0;
    }

  size_t w0_idx = 6 * (* i_pt), w1_idx = 6 * (* i_pt + 1);

  // if no Y-crossing, continue
  if(pt[w1_idx] * pt[w0_idx] > 0) return 0;

  int is_positive_vx = pt[w1_idx + 3] > 0.;
  if((ws -> use_positive_vx)^(is_positive_vx)) return 0;

  // Now record this point, using linear interpolation
  double
    y0  = pt[w0_idx + 1], y1  = pt[w1_idx + 1],
    py0 = pt[w0_idx + 4], py1 = pt[w1_idx + 4],
    z0  = pt[w0_idx + 2], z1  = pt[w1_idx + 2],
    pz0 = pt[w0_idx + 5], pz1 = pt[w1_idx + 5];

  double d = -pt[w0_idx] / (pt[w1_idx] - pt[w0_idx]);

  double
    ys  = (1. - d) * y0  + d * y1,
    pys = (1. - d) * py0 + d * py1,
    zs  = (1. - d) * z0  + d * z1,
    pzs = (1. - d) * pz0 + d * pz1;

  // DEBUG
  //printf("%f %f %f %f\n", ys, pys, zs, pzs);

  (ws -> points)[4 * (ws -> I_point)    ] = ys,
  (ws -> points)[4 * (ws -> I_point) + 1] = pys,
  (ws -> points)[4 * (ws -> I_point) + 2] = zs,
  (ws -> points)[4 * (ws -> I_point) + 3] = pzs;

  ++ (ws -> I_point);

  return 0;
}

int
_orbita_analyzer_poincare_section_x_kill(struct orbita_analyzer * wat)
{
  // free the workspace
  struct _orbita_analyzer_poincare_section_wspac * ws = 
    (struct _orbita_analyzer_poincare_section_wspac *)(wat -> wspac);
  free(ws -> points), free(ws);
  return 0;
}

struct orbita_analyzer *
orbita_analyzer_poincare_section_x(int N_pts, int N_max_cycles,
                                   int positive_vx)
{
  struct orbita_analyzer * wat = 
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  struct _orbita_analyzer_poincare_section_wspac * ws =
    (struct _orbita_analyzer_poincare_section_wspac *)
    malloc(sizeof(struct _orbita_analyzer_poincare_section_wspac));
  wat -> wspac  = (void *) ws;

  ws  -> points = (double *) malloc(sizeof(double) * 4 * N_pts),
  ws  -> N_points = N_pts,
  ws  -> N_cycles = N_max_cycles,
  ws  -> I_cycle  = 0,
  ws  -> I_point  = 0,
  ws  -> use_positive_vx = positive_vx;

  wat -> is_group = 0,
  wat -> N_analyzers = 0,
  wat -> group = NULL;

  wat -> peek = & _orbita_analyzer_poincare_section_x_peek;
  wat -> kill = & _orbita_analyzer_poincare_section_x_kill;

  return wat;
}

int
orbita_analyzer_poincare_section_dump(struct orbita_analyzer * wat,
                                      const char * filename)
{
  struct _orbita_analyzer_poincare_section_wspac * ws =
    (struct _orbita_analyzer_poincare_section_wspac *)(wat -> wspac);

  int N_pts = ws -> I_point;
  double * pts = ws -> points;

  // open a file and write the data block
  FILE * fp = fopen(filename, "wb");
  if(fp == NULL) return -1;
  fwrite((void *)pts, sizeof(double), N_pts * 4, fp);
  fclose(fp);

  return 0; 
}
