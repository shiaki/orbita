
#ifndef FERRERS_H

struct orbita_potential_ferrers2_paramset
{
  double a, b, c, rho, C;
  double W_i[20];
  /*
  double w000;
  double w100, w010, w001;
  double w110, w011, w101, w200, w020, w002;
  double w111, w120, w012, w201, w210, w021, w102;
  double w300, w030, w003;
  */
};

struct orbita_potential *
orbita_potential_ferrers2(double rho, double a, double b, double c);

#define FERRERS_H
#endif
