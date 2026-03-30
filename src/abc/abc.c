#include "abc.h"
#include "alg_inc.h"
#include "matrix/alg_matrix.h"
#include "memalloc/alg_memalloc.h"
#include "random/alg_random.h"
#include "vector/alg_vector.h"

abc_handle *abc_init(optim_handle optim, int pop_size, int max_count) {
    abc_handle *handle = ALG_MALLOC(sizeof(abc_handle));
    if (handle == NULL)
        return NULL;
    handle->population = alg_matrix_create(pop_size, optim.dim);
    handle->fitness = alg_vector_create(pop_size, 0.0);
    handle->count = ALG_CALLOC((size_t)pop_size, sizeof(int));
    if (handle->population == NULL || handle->fitness == NULL || handle->count == NULL) {
        alg_matrix_free(handle->population);
        alg_vector_free(handle->fitness);
        ALG_FREE(handle->count);
        ALG_FREE(handle);
        return NULL;
    }
    handle->optim = optim;
    handle->pop_size = pop_size;
    handle->max_count = max_count;
    alg_matrix_fill_random_vecs(handle->population, optim.l_range, optim.r_range, SET_ROW);
    optim_fresh(&handle->optim, handle->population, handle->fitness);
    return handle;
}

alg_state abc_fresh(abc_handle *handle, int gen) {
    alg_vector *cache_vec = alg_vector_create(handle->optim.dim, 0.0),
               *r_cache_vec = alg_vector_create(handle->optim.dim, 0.0),
               *copy_fitness = alg_vector_create_like(handle->fitness);
    for (int __iter = 0; __iter < gen; __iter++) {
        for (int pop = 0; pop < handle->pop_size; pop++) {
            alg_matrix_get_row(handle->population, cache_vec, pop);
            double cache_fitness = handle->optim.function(cache_vec);

            int index = alg_random_except_int(0, handle->pop_size, pop);
            alg_matrix_get_row(handle->population, r_cache_vec, index);
            for (int dim_iter = 0; dim_iter < handle->optim.dim; dim_iter++) {
                cache_vec->vector[dim_iter] +=
                    alg_random_float64(-1, 1) * (cache_vec->vector[dim_iter] - r_cache_vec->vector[dim_iter]);
            }
            alg_vector_claim_vecs(cache_vec, handle->optim.l_range, handle->optim.r_range);
            double r_cache_fitness = handle->optim.function(cache_vec);
            if (r_cache_fitness < cache_fitness) {
                alg_matrix_set_row(handle->population, pop, cache_vec);
                handle->count[pop] = 0;
            } else
                handle->count[pop] += 1;
        }

        optim_fresh(&handle->optim, handle->population, handle->fitness);

        double sum = alg_vector_sum(handle->fitness);
        copy_fitness->vector[0] = handle->fitness->vector[0] / sum;
        for (int pop = 1; pop < handle->pop_size; pop++) {
            copy_fitness->vector[pop] = handle->fitness->vector[pop] / sum + copy_fitness->vector[pop - 1];
        }
        double select_random = alg_random_float64(0, 1);
        int selected_index = handle->pop_size - 1;
        for (int pop = 0; pop < handle->pop_size; pop++) {
            if (select_random > copy_fitness->vector[pop])
                continue;
            selected_index = pop;
            break;
        }
        optim_fresh_best_solution(&handle->optim, handle->population, handle->fitness, selected_index);

        for (int pop = 0; pop < handle->pop_size; pop++) {
            if (handle->count[pop] < handle->max_count)
                continue;
            alg_matrix_get_row(handle->population, cache_vec, pop);
            for (int dim_iter = 0; dim_iter < handle->optim.dim; dim_iter++) {
                cache_vec->vector[dim_iter] = alg_random_float64(handle->optim.l_range->vector[dim_iter],
                                                                 handle->optim.r_range->vector[dim_iter]);
            }
            handle->count[pop] = 0;
        }
    }
    alg_vector_free(cache_vec);
    alg_vector_free(r_cache_vec);
    alg_vector_free(copy_fitness);
    return ALG_OK;
}

alg_state abc_free(abc_handle *handle) {
    alg_vector_free(handle->fitness);
    alg_matrix_free(handle->population);
    ALG_FREE(handle->count);
    ALG_FREE(handle);
    return ALG_OK;
}
