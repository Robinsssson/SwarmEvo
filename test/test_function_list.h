#ifndef __TEST_FUNCTION_LIST__H__
#define __TEST_FUNCTION_LIST__H__

#include "../src/optimization.h"

int test_abc(void);
int test_aco(void);
int test_aia(void);
int test_cs(void);
int test_cso(void);
int test_de(void);
int test_eda(void);
int test_ga(void);
int test_pso(void);
int test_sa(void);

#define TEST_FUNC 0
#define TEST_MAX_ITER 300
#define TEST_RET_VAL 2

#define VAL(v) v
#define VS 100.
#define DIM 6

#define TEST_OPTIM_DEFINE                                                                                              \
    optim_handle optim;                                                                                                \
    double *l_range = (double *)ALG_MALLOC(DIM * sizeof(double));                                                      \
    double *r_range = (double *)ALG_MALLOC(DIM * sizeof(double));                                                      \
    for (int i = 0; i < DIM; i++) {                                                                                    \
        l_range[i] = -VAL(VS);                                                                                         \
        r_range[i] = VAL(VS);                                                                                          \
    }                                                                                                                  \
    optim_init(&optim, DIM, test_function_list[TEST_FUNC], l_range, r_range);

#define TEST_OPTIM_END                                                                                                 \
    ALG_FREE(l_range);                                                                                                 \
    ALG_FREE(r_range);

#define CONCATENATE(x, y) x##y
#define EXPAND_AND_CONCATENATE(x, y) CONCATENATE(x, y)
#define HANDLE(method) EXPAND_AND_CONCATENATE(method, _handle)
#define INIT(method) EXPAND_AND_CONCATENATE(method, _init)
#define FRESH(method) EXPAND_AND_CONCATENATE(method, _fresh)
#define FREE(method) EXPAND_AND_CONCATENATE(method, _free)
#define TEST_FUNCTION(method) EXPAND_AND_CONCATENATE(test_, method)

#define TEST_FUNCTION_BODY(gen, iter_num)                                                                              \
    int TEST_FUNCTION(OPTIM_METHOD)(void) {                                                                            \
        TEST_OPTIM_DEFINE;                                                                                             \
        HANDLE(OPTIM_METHOD) *handle = INIT(OPTIM_METHOD)(optim, OPTIM_METHOD_ARGS);                                   \
        int i = 0;                                                                                                     \
        optim_print(&handle->optim);                                                                                   \
        do {                                                                                                           \
            FRESH(OPTIM_METHOD)(handle, gen);                                                                          \
            if (i % iter_num == 0)                                                                                     \
                optim_print(&handle->optim);                                                                           \
        } while (!FLOAT_COMPARE_IS(handle->optim.best_value, 0) && i++ < TEST_MAX_ITER);                               \
        optim_print(&handle->optim);                                                                                   \
        int converged = FLOAT_COMPARE_IS(handle->optim.best_value, TEST_RET_VAL);                                      \
        FREE(OPTIM_METHOD)(handle);                                                                                    \
        optim_free(&optim);                                                                                            \
        TEST_OPTIM_END;                                                                                                \
        int leaked = check_memory_leaks();                                                                             \
        return (converged && leaked == 0) ? TEST_PASSED : TEST_FAILED;                                                 \
    }

extern optimization test_function_list[];
extern double test_return_list[];
#endif //!__TEST_FUNCTION_LIST__H__
