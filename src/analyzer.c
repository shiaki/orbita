
#include "stdarg.h"
#include "stdlib.h"
#include "analyzer.h"

struct orbita_analyzer *
orbita_analyzer_group(int N, ...)
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  wat -> is_group    = 1;
  wat -> N_analyzers = N;
  wat -> group       = (struct orbita_analyzer **)
                           malloc(sizeof(struct orbita_analyzer *) * N);

  va_list  arglist;
  va_start(arglist, N);

  int i_wat;
  for(i_wat = 0; i_wat < N; ++ i_wat)
    (wat -> group)[i_wat] = va_arg(arglist, struct orbita_analyzer *);

  va_end(arglist);

  return wat;
}

int
orbita_analyzer_free(struct orbita_analyzer * wat)
{
  if(wat -> is_group)
    {
      int i_wat;
      for(i_wat = 0; i_wat < wat -> N_analyzers; ++ i_wat)
        orbita_analyzer_free((wat -> group)[i_wat]);
    }

  else (wat -> kill)(wat);

  free(wat);

  return 0;
}

int
orbita_analyzer_peek(struct orbita_analyzer * wat,
                     double * orb, int * I_tp, const int N)
{
  int st = 0; // status of integration

  if(wat -> is_group)
    {
      int i_wat;
      for(i_wat = 0; i_wat < wat -> N_analyzers; ++ i_wat)
        st = st | orbita_analyzer_peek((wat -> group)[i_wat], orb, I_tp, N);
    }

  else st = (wat -> peek)(wat -> wspac, orb, I_tp, N);

  return st;
}

int
_orbita_analyzer_dummy_peek(void * ws, double * pt,
                            int * i_pt, const int N_pts)
{
  if ((* i_pt) < 100) return 0; // <- that's for testing
  else return 0;
}

int
_orbita_analyzer_dummy_kill(struct orbita_analyzer * meow)
{
  return 0; // <- that's for testing
}

struct orbita_analyzer *
orbita_analyzer_dummy()
{
  struct orbita_analyzer * wat =
    (struct orbita_analyzer *) malloc(sizeof(struct orbita_analyzer));

  wat -> is_group    = 0,
  wat -> N_analyzers = 0,
  wat -> group       = NULL;

  wat -> peek  = & _orbita_analyzer_dummy_peek,
  wat -> wspac = NULL;

  wat -> kill  = & _orbita_analyzer_dummy_kill;

  return wat;
}
