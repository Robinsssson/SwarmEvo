#include "cso.h"
#include "alg_inc.h"
#include "math/alg_math.h"
#include "matrix/alg_matrix.h"
#include "memalloc/alg_memalloc.h"
#include "random/alg_random.h"
#include "vector/alg_vector.h"
#include <math.h>

cso_handle *cso_init(optim_handle optim, int pop_size, int rn, int hn, int cn, double fl) {
    if (pop_size != rn + hn + cn)
        return NULL;
    cso_handle *handle = ALG_MALLOC(sizeof(cso_handle));
    if (handle == NULL)
        return NULL;
    handle->fitness = alg_vector_create(pop_size, INFINITY);
    handle->position = alg_matrix_create(pop_size, optim.dim);
    handle->sorted_index = ALG_MALLOC(sizeof(int) * (size_t)pop_size);
    if (handle->fitness == NULL || handle->position == NULL || handle->sorted_index == NULL) {
        alg_vector_free(handle->fitness);
        alg_matrix_free(handle->position);
        ALG_FREE(handle->sorted_index);
        ALG_FREE(handle);
        return NULL;
    }
    handle->optim = optim;
    handle->pop_size = pop_size;
    handle->hn = hn; // 母鸡
    handle->rn = rn; // 公鸡
    handle->cn = cn;
    handle->fl = MATH_CLAIM(fl, 0, 2);

    alg_matrix_fill_random_vecs(handle->position, optim.l_range, optim.r_range, SET_ROW);
    optim_fresh(&handle->optim, handle->position, handle->fitness);
    alg_vector_sort_inplace(handle->fitness, handle->sorted_index, ALG_FALSE, alg_utils_greater);

    return handle;
}

static int rand_except(int left, int right, int except) {
    int index;
    do {
        index = alg_random_int(left, right);
    } while (except == index);
    return index;
}

static int find_parent(cso_handle *handle, int index) {
    int parent;
    if (index < handle->rn || index >= handle->pop_size)
        return -1;
    if (index >= handle->rn && index < handle->rn + handle->hn) {
        parent = index % handle->rn;
    } else {
        index -= handle->rn + handle->hn;
        parent = index % handle->hn + handle->rn;
    }
    return parent;
}

static int find_child(cso_handle *handle, int parent) {
    int index;
    do {
        index = handle->rn * alg_random_int(1, handle->hn / handle->rn + 2) + parent;
    } while (index >= handle->rn + handle->hn);
    return index;
}

alg_state cso_fresh(cso_handle *handle, int gen) {
    alg_vector *cache_vector = alg_vector_create(handle->optim.dim, 0.0),
               *k_vector = alg_vector_create(handle->optim.dim, 0.0),
               *r_vector = alg_vector_create(handle->optim.dim, 0.0);
    int index_k, iter_dim, group, i;
    double cache_fitness, k_fitness, r_fitness;
    double sigma, rand_val, c1, c2;
    for (int __iter = 0; __iter < gen; __iter++) {

        for (i = 0; i < handle->rn; i++) {
            alg_matrix_get_row(handle->position, cache_vector, handle->sorted_index[i]);
            cache_fitness = handle->optim.function(cache_vector);

            index_k = rand_except(0, handle->hn, handle->sorted_index[i]);
            alg_matrix_get_row(handle->position, k_vector, handle->sorted_index[index_k]);
            k_fitness = handle->optim.function(k_vector);
            sigma = k_fitness < cache_fitness
                        ? exp(alg_math_safe_divide(k_fitness - cache_fitness, fabs(cache_fitness) + 1e-9))
                        : 1;

            rand_val = alg_random_normal(0, sigma);
            for (iter_dim = 0; iter_dim < handle->optim.dim; iter_dim++)
                *alg_matrix_get_pos_mutval(handle->position, handle->sorted_index[i], iter_dim) *= 1 + rand_val;
            alg_matrix_clamp_vecs(handle->position, handle->optim.l_range, handle->optim.r_range, SET_ROW);
        }

        for (i = handle->rn; i < handle->rn + handle->hn; i++) {
            group = find_parent(handle, i);

            alg_matrix_get_row(handle->position, cache_vector, handle->sorted_index[i]);
            cache_fitness = handle->optim.function(cache_vector);

            alg_matrix_get_row(handle->position, k_vector, handle->sorted_index[group]);
            k_fitness = handle->optim.function(k_vector);

            alg_matrix_get_row(handle->position, r_vector,
                               handle->sorted_index[rand_except(0, handle->rn + handle->hn, i)]);
            r_fitness = handle->optim.function(r_vector);

            c1 = exp(alg_math_safe_divide(cache_fitness - k_fitness, fabs(cache_fitness) + 1e-9));
            c2 = exp(r_fitness - cache_fitness);
            if (isnan((float)c1) || isnan((float)c2) || isinf((float)c1) || isinf((float)c2)) {
                c1 = 0.0;
                c2 = 0.0;
            }
            for (iter_dim = 0; iter_dim < handle->optim.dim; iter_dim++)
                *alg_matrix_get_pos_mutval(handle->position, handle->sorted_index[i], iter_dim) +=
                    c1 * alg_random_float64(0, 1) * (k_vector->vector[iter_dim] - cache_vector->vector[iter_dim])
                    + c2 * alg_random_float64(0, 1) * (r_vector->vector[iter_dim] - cache_vector->vector[iter_dim]);

            alg_matrix_clamp_vecs(handle->position, handle->optim.l_range, handle->optim.r_range, SET_ROW);
        }
        for (i = handle->rn + handle->hn; i < handle->pop_size; i++) {
            index_k = find_child(handle, find_parent(handle, find_parent(handle, i)));
            for (iter_dim = 0; iter_dim < handle->optim.dim; iter_dim++)
                *alg_matrix_get_pos_mutval(handle->position, handle->sorted_index[i], iter_dim) +=
                    handle->fl
                    * (*alg_matrix_get_pos_val(handle->position, handle->sorted_index[index_k], iter_dim)
                       - *alg_matrix_get_pos_val(handle->position, handle->sorted_index[i], iter_dim));
        }
        alg_matrix_clamp_vecs(handle->position, handle->optim.l_range, handle->optim.r_range, SET_ROW);
        optim_fresh(&handle->optim, handle->position, handle->fitness);
        alg_vector_sort_inplace(handle->fitness, handle->sorted_index, ALG_FALSE, alg_utils_greater);
    }
    alg_vector_free(cache_vector);
    alg_vector_free(k_vector);
    alg_vector_free(r_vector);
    return ALG_OK;
}

alg_state cso_free(cso_handle *handle) {
    alg_matrix_free(handle->position);
    alg_vector_free(handle->fitness);
    ALG_FREE(handle->sorted_index);
    ALG_FREE(handle);
    return ALG_OK;
}
