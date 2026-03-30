#include "ga.h"
#include "alg_inc.h"
#include "matrix/alg_matrix.h"
#include "vector/alg_vector.h"
#include <math.h>
#include <stdio.h>

// 初始化遗传算法结构体
ga_handle *ga_init(optim_handle optim, int pop_size, double mutation_rate, double crossover_rate) {
    ga_handle *handle = ALG_MALLOC(sizeof(ga_handle));
    if (handle == NULL)
        return NULL;
    handle->optim = optim;
    handle->crossover_rate = crossover_rate; // 交叉率
    handle->mutation_rate = mutation_rate;   // 变异率
    handle->pop_size = pop_size;             // 种群大小

    // 初始化种群矩阵
    handle->population = alg_matrix_create(pop_size, optim.dim);
    if (handle->population == NULL) {
        ALG_FREE(handle);
        return NULL;
    }

    alg_matrix_fill_random_vecs(handle->population, optim.l_range, optim.r_range, SET_ROW); // 用随机值填充种群

    // 初始化适应度向量
    handle->fitness = alg_vector_create(handle->pop_size, 0.0);
    if (handle->fitness == NULL) {
        alg_matrix_free(handle->population);
        ALG_FREE(handle);
        return NULL;
    }

    optim_fresh(&handle->optim, handle->population, handle->fitness);
    return handle;
}

static alg_state sort_base_on_fitness(ga_handle *handle) {
    alg_vector *tmp_out = alg_vector_create_like(handle->fitness);
    if (tmp_out == NULL) {
        ERROR("THE VAL 'tmp_out' INIT ERROR");
        return ALG_ERROR;
    }
    int array[handle->pop_size];
    // base on fitness to copy in tmp_out and get index of array
    alg_vector_sort_copy(tmp_out, handle->fitness, array, alg_utils_greater);
    alg_matrix *copy_population = alg_matrix_copy(handle->population);
    if (copy_population == NULL) {
        ERROR("COPY POPULATION IS FAIL");
        alg_vector_free(tmp_out);
        return ALG_ERROR;
    }
    for (int i = 0; i < handle->pop_size; i++) {
        for (int j = 0; j < handle->optim.dim; j++) {
            alg_matrix_set_val(handle->population, i, j, *alg_matrix_get_pos_val(copy_population, array[i], j));
        }
    }
    for (int i = 0; i < tmp_out->size; i++)
        alg_vector_set_val(handle->fitness, i, tmp_out->vector[i]);
    alg_vector_free(tmp_out);
    alg_matrix_free(copy_population);
    optim_fresh(&handle->optim, handle->population, handle->fitness);
    return ALG_OK;
}

// 交叉操作：生成两个子代
static alg_state crossover(ga_handle *handle, const alg_vector *parent1, const alg_vector *parent2, alg_vector **child1,
                           alg_vector **child2) {
    int cross_number = alg_random_int(1, handle->optim.dim); // 随机选择交叉点
    alg_state state;
    // 将父代划分为两部分
    alg_vector *tmp_parent1_end = alg_vector_slice(parent1, cross_number, ALG_ALL_RANGE);
    alg_vector *tmp_parent2_end = alg_vector_slice(parent2, cross_number, ALG_ALL_RANGE);
    alg_vector *tmp_parent1_begin = alg_vector_slice(parent1, ALG_ALL_RANGE, cross_number);
    alg_vector *tmp_parent2_begin = alg_vector_slice(parent2, ALG_ALL_RANGE, cross_number);
    if (tmp_parent1_begin == NULL || tmp_parent2_begin == NULL) {
        ERROR("ERROR");
        return ALG_ERROR;
    }
    if (tmp_parent1_end == NULL || tmp_parent2_end == NULL || tmp_parent1_begin == NULL || tmp_parent2_begin == NULL) {
        ERROR("CREARE VECTOR SLICE IS ERROR");
        return ALG_ERROR;
    }

    // 交换部分基因，生成子代
    state = alg_vector_concat_inplace(&tmp_parent1_begin, tmp_parent2_end, ALG_VECTOR_CONCAT_R);
    if (state == ALG_ERROR) {
        ERROR("CHILD CREATE ERROR");
        return ALG_ERROR;
    }
    state = alg_vector_concat_inplace(&tmp_parent2_begin, tmp_parent1_end, ALG_VECTOR_CONCAT_R);
    if (state == ALG_ERROR) {
        ERROR("CHILD CREATE ERROR");
        return ALG_ERROR;
    }
    *child1 = tmp_parent1_begin;
    *child2 = tmp_parent2_begin;

    // 释放临时向量
    alg_vector_free(tmp_parent1_end);
    alg_vector_free(tmp_parent2_end);

    return ALG_OK;
}

// 变异操作：对子代进行变异
static void mutate(ga_handle *handle, alg_vector *vector) {
    if (alg_random_float64(0, 1) < handle->mutation_rate) {
        int mutation_point = alg_random_int(0, handle->optim.dim); // 随机选择变异点
        // 在变异点处用一个随机值进行变异
        alg_vector_set_val(vector, mutation_point,
                           alg_random_float64(handle->optim.l_range->vector[mutation_point],
                                              handle->optim.r_range->vector[mutation_point]));
    }
}

// 生成新一代种群
static alg_state generate_new_population(ga_handle *handle) {
    // 按照适应度对种群进行排序，选择优秀的父代
    sort_base_on_fitness(handle);

    // 选择父代的数量，这里我们选择偶数个父代
    int number_parents = (int)(round(handle->crossover_rate * handle->pop_size) / 2) * 2;
    // 新种群的矩阵，大小是父代数量的两倍（每对父代生成两个子代）
    alg_matrix *new_population = alg_matrix_create(number_parents * 2, handle->optim.dim);
    if (new_population == NULL) {
        ERROR("NEW POPULATION CREATE ERROR");
        return ALG_ERROR;
    }

    // 用于存放子代的数组
    alg_vector **child_list = (alg_vector **)ALG_MALLOC((size_t)(2 * number_parents) * sizeof(alg_vector *));
    if (child_list == NULL) {
        ERROR("CHILD LIST ALLOC ERROR");
        alg_matrix_free(new_population);
        return ALG_ERROR;
    }
    for (int i = 0; i < 2 * number_parents; i++) {
        child_list[i] = NULL;
    }
    // 开始交叉并生成子代
    for (int i = 0; i < number_parents; i += 2) {
        // 从种群中选择两个父代
        alg_vector *parent1 = alg_vector_from_matrix_row(handle->population, i);
        alg_vector *parent2 = alg_vector_from_matrix_row(handle->population, i + 1);
        if (parent1 == NULL || parent2 == NULL) {
            alg_vector_free(parent1);
            alg_vector_free(parent2);
            for (int k = 0; k < i; k++)
                alg_vector_free(child_list[k]);
            ALG_FREE(child_list);
            alg_matrix_free(new_population);
            return ALG_ERROR;
        }
        // 进行交叉生成两个子代
        alg_vector *child1 = NULL, *child2 = NULL;
        if (crossover(handle, parent1, parent2, &child1, &child2) != ALG_OK) {
            alg_vector_free(parent1);
            alg_vector_free(parent2);
            for (int k = 0; k < i; k++)
                alg_vector_free(child_list[k]);
            ALG_FREE(child_list);
            alg_matrix_free(new_population);
            return ALG_ERROR;
        }

        // 进行变异操作
        mutate(handle, child1);
        mutate(handle, child2);

        // 存储子代
        child_list[i] = child1;
        child_list[i + 1] = child2;

        // 释放父代向量
        alg_vector_free(parent1);
        alg_vector_free(parent2);
    }

    // 将交叉和变异后的子代填充到新种群
    for (int i = 0; i < number_parents; i += 2) {
        for (int j = 0; j < handle->optim.dim; j++) {
            alg_matrix_set_val(new_population, i, j, child_list[i]->vector[j]);
            alg_matrix_set_val(new_population, i + 1, j, child_list[i + 1]->vector[j]);
        }
    }

    // 选择精英策略：将适应度最好的父代保存到新种群中
    for (int i = number_parents; i < 2 * number_parents; i++) {
        for (int j = 0; j < handle->optim.dim; j++) {
            // 保留精英父代（最好适应度的个体）
            alg_matrix_set_val(new_population, i, j,
                               *alg_matrix_get_pos_val(handle->population, i - number_parents, j));
        }
    }

    // 替换旧种群为新种群
    alg_matrix_free(handle->population);
    for (int i = 0; i < 2 * number_parents; i++) {
        alg_vector_free(child_list[i]);
    }
    ALG_FREE(child_list);
    handle->population = new_population;
    handle->pop_size = handle->population->row;
    return ALG_OK;
}

// 主遗传算法循环，更新种群
alg_state ga_fresh(ga_handle *handle, int gen) {
    for (int __iter = 0; __iter < gen; __iter++) {
        optim_fresh(&handle->optim, handle->population, handle->fitness);
        generate_new_population(handle); // 生成新一代种群
    }
    return ALG_OK;
}

// 释放遗传算法的内存
alg_state ga_free(ga_handle *ga) {
    if (ga == NULL)
        return ALG_ERROR;
    alg_matrix_free(ga->population); // 释放种群矩阵
    alg_vector_free(ga->fitness);    // 释放适应度向量
    ALG_FREE(ga);                    // 释放GA结构体内存
    return ALG_OK;
}
