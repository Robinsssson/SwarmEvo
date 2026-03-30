#include "aia.h"
#include "alg_inc.h"
#include "matrix/alg_matrix.h"
#include "vector/alg_vector.h"
#include <math.h>

aia_handle *aia_init(optim_handle optim, int num_antibodies, int clone_factor, double mutation_rate) {
    aia_handle *handle = ALG_MALLOC(sizeof(aia_handle));
    if (handle == NULL) {
        ERROR("ERROR INIT IN AIA_HANDLE OBJECT");
        return NULL;
    }
    handle->population = alg_matrix_create(num_antibodies, optim.dim);
    if (handle->population == NULL) {
        ERROR("ERROR INIT IN POPULATION OBJECT");
        ALG_FREE(handle);
        return NULL;
    }
    handle->fitness = alg_vector_create(num_antibodies, 0.0);
    if (handle->fitness == NULL) {
        ERROR("ERROR INIT IN FITNESS OBJECT");
        alg_matrix_free(handle->population);
        ALG_FREE(handle);
        return NULL;
    }
    alg_matrix_fill_random_vecs(handle->population, optim.l_range, optim.r_range, SET_ROW);
    handle->num_antibodies = num_antibodies;
    handle->clone_factor = clone_factor;
    handle->mutation_rate = mutation_rate;
    handle->optim = optim;
    optim_fresh(&handle->optim, handle->population, handle->fitness);
    return handle;
}

static alg_vector *evaluate_fitness(aia_handle *handle, alg_matrix *mat) {
    alg_vector *_tmp_vec = alg_vector_create(handle->optim.dim, 0.0);
    alg_vector *vecs = alg_vector_create(mat->row, 0.0);
    for (int row = 0; row < mat->row; row++) {
        alg_matrix_get_row(mat, _tmp_vec, row);
        alg_vector_set_val(vecs, row, handle->optim.function(_tmp_vec));
    }
    alg_vector_free(_tmp_vec);
    return vecs;
}

static alg_matrix *clone_and_mutate(aia_handle *handle, alg_vector *best_sol) {
    int len = handle->num_antibodies;
    int number_clones[len];

    double min_fit = handle->fitness->vector[0];
    double max_fit = handle->fitness->vector[0];
    for (int i = 1; i < handle->fitness->size; i++) {
        if (handle->fitness->vector[i] < min_fit)
            min_fit = handle->fitness->vector[i];
        if (handle->fitness->vector[i] > max_fit)
            max_fit = handle->fitness->vector[i];
    }
    double range_fit = max_fit - min_fit + 1e-6;

    double sum_val = 0;
    for (int i = 0; i < len; i++) {
        number_clones[i] = (int)ceil(handle->clone_factor * (1.0 - (handle->fitness->vector[i] - min_fit) / range_fit));
        if (number_clones[i] < 1)
            number_clones[i] = 1;
        sum_val += number_clones[i];
    }

    alg_matrix *clones = alg_matrix_create((int)sum_val, handle->optim.dim);
    alg_vector *clone = alg_vector_create(handle->optim.dim, 0.0);
    int count = 0;
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < number_clones[i]; j++) {
            alg_matrix_get_row(handle->population, clone, i);
            double mutation_strength = handle->mutation_rate * (1.0 - (double)j / number_clones[i]);
            for (int iter = 0; iter < clone->size; iter++) {
                double range_dim = handle->optim.r_range->vector[iter] - handle->optim.l_range->vector[iter];
                double r = alg_random_float64(-1, 1);
                clone->vector[iter] = best_sol->vector[iter] + range_dim * mutation_strength * r * r * (r > 0 ? 1 : -1);
            }
            alg_vector_claim_vecs(clone, handle->optim.l_range, handle->optim.r_range);
            alg_matrix_set_row(clones, count, clone);
            count++;
        }
    }
    alg_vector_free(clone);
    return clones;
}

static void select_new_population(aia_handle *handle, alg_matrix **clones) {
    alg_matrix_concat(clones, handle->population, CONCAT_AXIS_DY);
    alg_vector *fitness = evaluate_fitness(handle, *clones);
    alg_vector *copy_fitness = alg_vector_create_like(fitness);
    int sorted_indices[fitness->size];
    alg_vector_sort_copy(copy_fitness, fitness, sorted_indices, alg_utils_greater);
    for (int i = 0; i < handle->population->row; i++)
        for (int j = 0; j < handle->population->col; j++) {
            alg_matrix_set_val(handle->population, i, j, *alg_matrix_get_pos_val(*clones, sorted_indices[i], j));
        }
    alg_vector_free(fitness);
    alg_vector_free(copy_fitness);
}

alg_state aia_fresh(aia_handle *handle, int gen) {
    alg_vector *best_sol = alg_vector_create(handle->optim.dim, 0.0);
    double best_val = INFINITY;

    for (int __iter = 0; __iter < gen; __iter++) {
        alg_matrix *clones = clone_and_mutate(handle, best_sol);
        select_new_population(handle, &clones);
        alg_matrix_free(clones);
        optim_fresh(&handle->optim, handle->population, handle->fitness);

        for (int i = 0; i < handle->num_antibodies; i++) {
            if (handle->fitness->vector[i] < best_val) {
                best_val = handle->fitness->vector[i];
                for (int d = 0; d < handle->optim.dim; d++)
                    best_sol->vector[d] = *alg_matrix_get_pos_val(handle->population, i, d);
            }
        }

        double min_f = handle->fitness->vector[0], max_f = handle->fitness->vector[0];
        for (int i = 1; i < handle->fitness->size; i++) {
            if (handle->fitness->vector[i] < min_f)
                min_f = handle->fitness->vector[i];
            if (handle->fitness->vector[i] > max_f)
                max_f = handle->fitness->vector[i];
        }
        if (max_f - min_f < 5.0) {
            for (int idx = 0; idx < 3; idx++) {
                int random_idx = alg_random_int(0, handle->num_antibodies);
                for (int d = 0; d < handle->optim.dim; d++) {
                    double range = handle->optim.r_range->vector[d] - handle->optim.l_range->vector[d];
                    double val = best_sol->vector[d] + alg_random_float64(-0.1, 0.1) * range;
                    val = MATH_CLAIM(val, handle->optim.l_range->vector[d], handle->optim.r_range->vector[d]);
                    alg_matrix_set_val(handle->population, random_idx, d, val);
                }
            }
            optim_fresh(&handle->optim, handle->population, handle->fitness);
        }
    }

    // 最终局部搜索 - 坐标下降法
    alg_vector *tmp_sol = alg_vector_create(handle->optim.dim, 0.0);
    double step = 0.01;
    for (int round = 0; round < 5; round++) {
        for (int d = 0; d < handle->optim.dim; d++) {
            for (int dir = -1; dir <= 1; dir += 2) {
                for (int i = 0; i < 20; i++) {
                    for (int j = 0; j < handle->optim.dim; j++)
                        tmp_sol->vector[j] = best_sol->vector[j];
                    tmp_sol->vector[d] += dir * step;
                    alg_vector_claim_vecs(tmp_sol, handle->optim.l_range, handle->optim.r_range);
                    double new_val = handle->optim.function(tmp_sol);
                    if (new_val < best_val) {
                        best_val = new_val;
                        for (int k = 0; k < handle->optim.dim; k++)
                            best_sol->vector[k] = tmp_sol->vector[k];
                    } else {
                        break;
                    }
                }
            }
        }
        step *= 0.5;
    }
    for (int d = 0; d < handle->optim.dim; d++)
        handle->optim.best_solution->vector[d] = best_sol->vector[d];
    handle->optim.best_value = best_val;
    alg_vector_free(tmp_sol);
    alg_vector_free(best_sol);
    return ALG_OK;
}

alg_state aia_free(aia_handle *handle) {
    alg_vector_free(handle->fitness);
    alg_matrix_free(handle->population);
    ALG_FREE(handle);
    return ALG_OK;
}
