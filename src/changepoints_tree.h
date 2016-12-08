#ifndef CHANGEPOINTS_TREE_H
#define CHANGEPOINTS_TREE_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "contrasts.h"

#define CPTS_LEN_STEP 128

typedef struct cpt_tree_node{
  
  double max;
  unsigned int *index;
  unsigned int n_intervals;
  unsigned cpt;
  struct cpt_tree_node *left_node, *right_node;
  
} cpt_tree_node_t; 

typedef struct cpts{
  
  unsigned int *cpt; 
  unsigned int n_cpt;
  double min_max;
  
  cpt_tree_node_t *node;
  cpt_tree_node_t *parent_node;
  
} cpts_t;

typedef struct solution_path{
  
  cpts_t *cpts;
  double *th;
  unsigned int n_th;
  
  
} solution_path_t;


void build_tree(cpt_tree_node_t **node, cpt_tree_node_t **parent_node, unsigned int start, unsigned int end,
                double th, contrasts_t *contrasts, eval_contrast_fun_t eval_contrast_fun);
void get_changepoints(cpt_tree_node_t **node, cpts_t *cpts, unsigned int start, unsigned int end, unsigned int min_dist);
solution_path_t *solution_path(contrasts_t *contrasts, eval_contrast_fun_t eval_contrast_fun, unsigned int min_dist);

void destroy_tree(cpt_tree_node_t **node);
void destroy_solution_path(solution_path_t **solution_path);

int compare_cpts_t(const cpts_t *a, const cpts_t *b, unsigned int n_obs);

#endif 