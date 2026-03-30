#include "basic_opti.h"
#include "alg_inc.h"
#include "matrix/alg_matrix.h"
#include "memalloc/alg_memalloc.h"
#include "vector/alg_vector.h"
#include <math.h>
#include <stdio.h>

alg_state optim_init(optim_handle *handle, int dim, optimization function, double l_range[], double r_range[]) {
    handle->function = function;
    handle->l_range = alg_vector_create(dim, 0.0);
    handle->r_range = alg_vector_create(dim, 0.0);
    handle->best_solution = alg_vector_create(dim, INFINITY);
    for (int i = 0; i < dim; i++) {
        handle->l_range->vector[i] = l_range[i];
        handle->r_range->vector[i] = r_range[i];
    }
    handle->dim = dim;
    return ALG_OK;
}

alg_state optim_free(optim_handle *handle) {
    alg_vector_free(handle->l_range);
    alg_vector_free(handle->r_range);
    alg_vector_free(handle->best_solution);
    return ALG_OK;
}

alg_state optim_fresh(optim_handle *handle, alg_matrix *population, alg_vector *fitness) {
    alg_vector *tmp = alg_vector_create(handle->dim, 0.0);
    for (int i = 0; i < population->row; i++) {
        for (int j = 0; j < population->col; j++) {
            tmp->vector[j] = *alg_matrix_get_pos_val(population, i, j);
        }
        fitness->vector[i] = handle->function(tmp);
    }
    int index = alg_vector_compare_val(fitness, alg_utils_greater);
    optim_fresh_best_solution(handle, population, fitness, index);
    alg_vector_free(tmp);
    return ALG_OK;
}

alg_state optim_fresh_best_solution(optim_handle *handle, alg_matrix *population, alg_vector *fitness, int index) {
    handle->best_value = fitness->vector[index];
    for (int i = 0; i < population->col; i++)
        handle->best_solution->vector[i] = *alg_matrix_get_pos_val(population, index, i);
    return ALG_OK;
}

alg_state optim_print(optim_handle *handle) {
    char *str = alg_vector_print_str(handle->best_solution);
    printf("OPTIMIZATION best value %.4f, %s\n", handle->best_value, str);
    ALG_FREE(str);
    return ALG_OK;
}
