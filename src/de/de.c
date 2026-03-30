#include "de.h"
#include "alg_inc.h"
#include "matrix/alg_matrix.h"
#include "memalloc/alg_memalloc.h"
#include "random/alg_random.h"
#include "vector/alg_vector.h"
#include <math.h>
#include <stdio.h>

static alg_boolean check_if(int *idx, int len, int index) {
    for (int i = 0; i < len; i++)
        if (idx[i] == index)
            return ALG_TRUE;
    return ALG_FALSE;
}

static alg_state get_from_index(const alg_matrix *mat, double *arr, int length, int index) {
    if (length == mat->col) {
        for (int i = 0; i < length; i++) {
            arr[i] = *alg_matrix_get_pos_val(mat, index, i);
        }
    } else if (length == mat->row) {
        for (int i = 0; i < length; i++) {
            arr[i] = *alg_matrix_get_pos_val(mat, i, index);
        }
    } else
        return ALG_ERROR;
    return ALG_OK;
}

de_handle *de_init(optim_handle optim, int pop_size, double f, double cr) {
    de_handle *handle = ALG_MALLOC(sizeof(de_handle));
    if (handle == NULL)
        return NULL;
    handle->fitness = alg_vector_create(pop_size, INFINITY);
    handle->population = alg_matrix_create(pop_size, optim.dim);
    if (handle->fitness == NULL || handle->population == NULL) {
        alg_vector_free(handle->fitness);
        alg_matrix_free(handle->population);
        ALG_FREE(handle);
        return NULL;
    }
    handle->optim = optim;
    alg_matrix_fill_random_vecs(handle->population, handle->optim.l_range, handle->optim.r_range, SET_ROW);
    handle->pop_size = pop_size;
    handle->f = f;
    handle->cr = cr;
    optim_fresh(&handle->optim, handle->population, handle->fitness);
    return handle;
}

alg_state de_fresh(de_handle *handle, int gen) {
    int idx[3];
    int j_rand;
    alg_vector *u = alg_vector_create(handle->optim.dim, 0.0);
    double u_fitness;
    double x1[handle->optim.dim], x2[handle->optim.dim], x3[handle->optim.dim], v[handle->optim.dim];
    for (int __iter = 0; __iter < gen; __iter++) {
        for (int i = 0; i < handle->pop_size; i++) {
            do {
                alg_random_sample_unique(0, (int)handle->pop_size, 3, idx);
            } while (check_if(idx, 3, i) == ALG_TRUE);

            if (get_from_index(handle->population, x1, handle->optim.dim, idx[0]) != ALG_OK
                || get_from_index(handle->population, x2, handle->optim.dim, idx[1]) != ALG_OK
                || get_from_index(handle->population, x3, handle->optim.dim, idx[2]) != ALG_OK) {
                return ALG_ERROR;
            }

            for (int i_dim = 0; i_dim < handle->optim.dim; i_dim++)
                v[i_dim] = x1[i_dim] + handle->f * (x2[i_dim] - x3[i_dim]);

            // 交叉操作
            alg_matrix_get_row(handle->population, u, i);
            j_rand = alg_random_int(0, handle->optim.dim);
            for (int j = 0; j < handle->optim.dim; j++) {
                if (alg_random_float64(0, 1) <= handle->cr || j == j_rand)
                    alg_vector_set_val(u, j, v[j]);
            }
            alg_vector_claim_vecs(u, handle->optim.l_range, handle->optim.r_range);
            // 选择操作
            u_fitness = handle->optim.function(u);
            if (u_fitness < handle->fitness->vector[i]) {
                alg_matrix_set_row(handle->population, i, u);
                handle->fitness->vector[i] = u_fitness;
            }
        }
        optim_fresh(&handle->optim, handle->population, handle->fitness);
    }
    alg_vector_free(u);
    return ALG_OK;
}

alg_state de_free(de_handle *handle) {
    alg_vector_free(handle->fitness);
    alg_matrix_free(handle->population);
    ALG_FREE(handle);
    return ALG_OK;
}
