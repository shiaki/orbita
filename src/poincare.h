
#ifndef POINCARE_H

struct _orbita_analyzer_poincare_section_wspac
{
  int N_points, N_cycles;
  int I_point,  I_cycle;
  int use_positive_vx;

  double * points;
};

int _orbita_analyzer_poincare_section_x_peek(void * wspac,
      double * pt, int * i_pt, const int N_pts);

int _orbita_analyzer_poincare_section_x_kill(struct orbita_analyzer * wat);

struct orbita_analyzer *
orbita_analyzer_poincare_section_x(int N_pts, int N_max_cycles,
                                   int positive_vx);

int orbita_analyzer_poincare_section_dump(struct orbita_analyzer * wat,
                                      const char * filename);
#define POINCARE_H
#endif
