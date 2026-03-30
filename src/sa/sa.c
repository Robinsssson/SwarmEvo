#include "sa.h"
#include "math/alg_math.h"
#include "vector/alg_vector.h"

static void sa_fresh_best_solution(sa_handle *handle) {
    for (int i = 0; i < handle->optim.dim; i++) {
        handle->optim.best_solution->vector[i] = handle->best_solution->vector[i];
    }
    handle->optim.best_value = handle->best_energy;
}

sa_handle *sa_init(optim_handle optim, double temperature, double cooling_rate) {
    sa_handle *handle = ALG_MALLOC(sizeof(sa_handle));
    if (handle == NULL) {
        ERROR("SA HANDLE INIT ERROR");
        return NULL;
    }
    handle->current_solution = alg_vector_create(optim.dim, 0.0);
    for (int i = 0; i < optim.dim; i++) {
        handle->current_solution->vector[i] = alg_random_float64(optim.l_range->vector[i], optim.r_range->vector[i]);
    }
    handle->best_solution = alg_vector_create_like(handle->current_solution);

    handle->new_solution = alg_vector_create_like(handle->current_solution);

    handle->optim = optim;
    handle->temperature = temperature;
    handle->cooling_rate = cooling_rate;
    handle->current_energy = optim.function(handle->current_solution);
    handle->best_energy = handle->current_energy;
    sa_fresh_best_solution(handle);
    return handle;
}

alg_state sa_fresh(sa_handle *handle, int gen) {
    for (int __iter = 0; __iter < gen; __iter++) {
        for (int i = 0; i < handle->optim.dim; i++) {
            handle->new_solution->vector[i] = handle->current_solution->vector[i] + alg_random_float64(-1, 1);
        }
        alg_vector_claim_vecs(handle->new_solution, handle->optim.l_range, handle->optim.r_range);

        double new_energy = handle->optim.function(handle->new_solution);

        if (new_energy < handle->current_energy
            || exp(alg_math_safe_divide(handle->current_energy - new_energy, handle->temperature + 1e-6))
                   > alg_random_float64(0, 1)) {
            handle->current_energy = new_energy;
            for (int i = 0; i < handle->optim.dim; i++)
                handle->current_solution->vector[i] = handle->new_solution->vector[i];
        }

        if (handle->current_energy < handle->best_energy) {
            handle->best_energy = handle->current_energy;
            for (int i = 0; i < handle->optim.dim; i++)
                handle->best_solution->vector[i] = handle->current_solution->vector[i];
        }
        handle->temperature *= handle->cooling_rate;
        handle->temperature = handle->temperature > 1 ? handle->temperature : 1;
        sa_fresh_best_solution(handle);
    }
    return ALG_OK;
}

alg_state sa_free(sa_handle *handle) {
    alg_vector_free(handle->best_solution);
    alg_vector_free(handle->current_solution);
    alg_vector_free(handle->new_solution);
    ALG_FREE(handle);
    return ALG_OK;
}
