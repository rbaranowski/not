#ifndef NOT_H
#define NOT_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "contrasts.h"
#include "changepoints_tree.h"

#define IDX(i,j,ld) ((((j)-1) * (ld))+((i)-1))

typedef struct ip_max{
  unsigned int arg_max;
  double max;
  double abs_max;
} ip_max_t;

typedef struct ips {
  unsigned int *index;
  unsigned int *s;
  unsigned int *e;
  unsigned int *cpt;
  double *max;
  double *abs_max;
  unsigned int M;
  unsigned int n;
} ips_t;


typedef struct notres{
  unsigned int *s;
  unsigned int *e;
  unsigned int *cpt;
  double *max;
  double *minth;
  unsigned int *scale;
  unsigned int M;
  unsigned int n;
} not_res_t;

SEXP not_r_wrapper(SEXP x, SEXP intervals, SEXP method, SEXP contrast_type, SEXP parallel, SEXP augmented);
SEXP solution_path_t_to_list(solution_path_t *solution_path);
SEXP contrasts_t_to_dataframe(contrasts_t *contrasts);


#endif