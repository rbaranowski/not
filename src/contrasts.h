#ifndef CONTRASTS_H
#define CONTRASTS_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

typedef struct max_contrast{
  unsigned int arg_max;
  double max;
} max_contrast_t;

typedef struct contrasts{
  unsigned int *index;
  unsigned int *start;
  unsigned int *end;
  unsigned int *length;
  unsigned int *arg_max;
  double *max;
  unsigned int n_intervals;
  double *x;
  unsigned int n_obs;
} contrasts_t;


typedef max_contrast_t (*eval_contrast_fun_t)(double *, unsigned int);

#define SQRT_THREE	1.73205080756888
#define SQRT_FIVE	  2.23606797749979
#define SQRT_SIX	  2.44948974278318
#define SQRT_TWELVE	  3.46410161513775

void alloc_contrasts(contrasts_t **contrasts, unsigned int n_intervals, double *x, unsigned int n_obs);
void destroy_contrasts(contrasts_t **contrasts);

max_contrast_t intercept_contrast(double *x, unsigned int n_obs);
max_contrast_t intercept_signs_contrast(double *x, unsigned int n_obs);
max_contrast_t intercept_and_slope_contrast(double *x, unsigned int n_obs);
max_contrast_t slope_contrast(double *x, unsigned int n_obs);
max_contrast_t intercept_slope_and_quadratic_contrast(double *x, unsigned int n_obs);
max_contrast_t intercept_and_volatility_contrast(double *x, unsigned int n_obs);


contrasts_t *eval_contrasts(double *x, unsigned int n_obs, unsigned int *intervals, unsigned int n_intervals,
                            eval_contrast_fun_t eval_contrast_fun, unsigned int parallel);



#endif