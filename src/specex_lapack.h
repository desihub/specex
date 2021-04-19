#ifndef SPECEX_LAPACK_H_
#define SPECEX_LAPACK_H_

extern "C" {
  int specex_posv(int, const double *, const double *);
  int specex_potri(int, const double *);  
}

#endif
