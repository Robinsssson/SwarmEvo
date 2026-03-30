#ifndef __BASIC_OPTI_H
#define __BASIC_OPTI_H
#include "../algmath/algmath.h"

/*
 * function pointer for optimization
 */
typedef double (*optimization)(alg_vector *);
typedef struct {
    int dim;
    alg_vector *l_range, *r_range;
    optimization function;
    alg_vector *best_solution;
    double best_value;
} optim_handle;

ALG_MATH_API alg_state optim_init(optim_handle *handle, int dim, optimization function, double l_range[],
                                  double r_range[]);
ALG_MATH_API alg_state optim_free(optim_handle *handle);
ALG_MATH_API alg_state optim_fresh(optim_handle *handle, alg_matrix *population, alg_vector *fitness);
ALG_MATH_API alg_state optim_fresh_best_solution(optim_handle *handle, alg_matrix *population, alg_vector *fitness,
                                                 int index);
ALG_MATH_API alg_state optim_print(optim_handle *handle);
#endif //!__BASIC_OPTI_H
