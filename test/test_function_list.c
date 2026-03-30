#include "test_function_list.h"
static double sphere_function(alg_vector *vec) {
    double ret = 0.0;
    for (int i = 0; i < vec->size; i++) {
        ret += pow(vec->vector[i], 2);
    }
    return ret + 2;
}

static double rosenbrock(alg_vector *vec) {
    double ret = 0.0;
    for (int i = 0; i < vec->size - 1; i++) {
        ret += pow(1 - vec->vector[i], 2) + 100 * pow(vec->vector[i + 1] - pow(vec->vector[i], 2), 2);
    }
    return ret;
}

optimization test_function_list[] = {
    sphere_function,
    rosenbrock,
};

double test_return_list[] = {
    2.0, // sphere_function minimum is 0 + 2 = 2
    0.0, // rosenbrock minimum is 0
};
