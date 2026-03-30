#include "pso.h"
#include "alg_inc.h"
#include "algmath.h"
#include "matrix/alg_matrix.h"
#include "vector/alg_vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Helper function to ALG_FREE PSO resources
static void pso_free_internal(pso_handle *handle) {
    alg_matrix_free(handle->g_best);
    alg_matrix_free(handle->p_best);
    alg_matrix_free(handle->position);
    alg_matrix_free(handle->vec);
    alg_vector_free(handle->fitness);
    ALG_FREE(handle->p_best_fitness);
}

pso_handle *pso_init(optim_handle optim, int pop_size, double w, double c1, double c2) {
    pso_handle *handle = ALG_MALLOC(sizeof(pso_handle));
    if (!handle)
        return NULL;

    handle->vec = alg_matrix_create(pop_size, optim.dim);
    handle->position = alg_matrix_create(pop_size, optim.dim);
    handle->g_best = alg_matrix_create(1, optim.dim);
    handle->fitness = alg_vector_create(pop_size, 0.0);
    if (!handle->vec || !handle->position || !handle->g_best || !handle->fitness) {
        pso_free_internal(handle);
        return NULL;
    }

    alg_matrix_fill_random_vecs(handle->position, optim.l_range, optim.r_range, SET_ROW);
    handle->p_best = alg_matrix_copy(handle->position);
    if (!handle->p_best) {
        pso_free_internal(handle);
        return NULL;
    }
    alg_matrix_fill_random(handle->vec, -1, 1);

    handle->p_best_fitness = ALG_MALLOC((size_t)pop_size * sizeof(double));
    if (!handle->p_best_fitness) {
        pso_free_internal(handle);
        return NULL;
    }
    for (int i = 0; i < pop_size; i++)
        handle->p_best_fitness[i] = INFINITY;

    handle->optim = optim;
    handle->g_best_index = 0;
    handle->g_best_fitness = INFINITY;
    handle->c1 = c1;
    handle->c2 = c2;
    handle->pop_size = pop_size;
    handle->w = w;
    optim_fresh(&handle->optim, handle->position, handle->fitness);
    return handle;
}

static double pso_calculate_fitness(pso_handle *handle, int index) {
    alg_vector *vec = alg_vector_from_matrix_row(handle->position, index);
    double fitness = handle->optim.function(vec);
    alg_vector_free(vec);
    return fitness;
}

alg_state pso_fresh(pso_handle *handle, int gen) {
    alg_vector *current_position = alg_vector_create(handle->optim.dim, 0.0);
    alg_vector *best_position = alg_vector_create(handle->optim.dim, 0.0);
    for (int __iter = 0; __iter < gen; __iter++) {
        double best_val = INFINITY;
        int best_idx = 0;

        for (int i = 0; i < handle->pop_size; i++) {
            double fitness = pso_calculate_fitness(handle, i);

            if (fitness < handle->p_best_fitness[i]) {
                handle->p_best_fitness[i] = fitness;
                alg_matrix_get_row(handle->position, current_position, i);
                alg_matrix_set_row(handle->p_best, i, current_position);
            }

            if (fitness < best_val) {
                best_val = fitness;
                best_idx = i;
            }
        }

        if (best_val < handle->g_best_fitness) {
            alg_matrix_get_row(handle->position, best_position, best_idx);
            alg_matrix_set_row(handle->g_best, 0, best_position);
            handle->g_best_fitness = best_val;
            handle->g_best_index = best_idx;
        }

        alg_matrix_get_row(handle->g_best, best_position, 0);

        for (int i = 0; i < handle->pop_size; i++) {
            double r1 = alg_random_float64(0, 1);
            double r2 = alg_random_float64(0, 1);

            for (int d = 0; d < handle->optim.dim; d++) {
                double *v = alg_matrix_get_pos_mutval(handle->vec, i, d);
                double pos = *alg_matrix_get_pos_val(handle->position, i, d);
                double pbest = *alg_matrix_get_pos_val(handle->p_best, i, d);

                *v = handle->w * (*v) + handle->c1 * r1 * (pbest - pos)
                     + handle->c2 * r2 * (best_position->vector[d] - pos);
            }
        }

        double max_vel = 0.2;
        alg_matrix_clamp(handle->vec, -max_vel, max_vel);

        alg_matrix_add_inplace(handle->position, handle->vec);
        alg_matrix_clamp_vecs(handle->position, handle->optim.l_range, handle->optim.r_range, SET_ROW);

        optim_fresh(&handle->optim, handle->position, handle->fitness);
    }
    alg_vector_free(best_position);
    alg_vector_free(current_position);
    return ALG_OK;
}

alg_state pso_free(pso_handle *handle) {
    pso_free_internal(handle);
    ALG_FREE(handle);
    return ALG_OK;
}
