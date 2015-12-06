
#ifndef ANALYZER_H

struct orbita_analyzer
{
  int    is_group;
  int    N_analyzers;
  struct orbita_analyzer ** group;

  int (* peek)(void *, double *, int *, const int);
  void * wspac;

  int (* kill)(struct orbita_analyzer *);
};

struct orbita_analyzer * orbita_analyzer_group(int N, ...);

int orbita_analyzer_free(struct orbita_analyzer *);

int orbita_analyzer_peek(struct orbita_analyzer *,
    double *, int *, const int);

struct orbita_analyzer * orbita_analyzer_dummy();

#define ANALYZER_H
#endif
