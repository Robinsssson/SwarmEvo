#ifndef DEBUG_MEMORY_H
#define DEBUG_MEMORY_H

#include <stddef.h> // 为 size_t 提供定义

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 分配内存并记录信息，用于检测内存泄漏。
 *
 * @param size 要分配的内存大小。
 * @return 分配的内存块指针，失败返回 NULL。
 */
void *debug_malloc(size_t size);

/**
 * @brief 分配并清零内存块，用于检测内存泄漏。
 *
 * @param num 要分配的元素个数。
 * @param size 每个元素的大小。
 * @return 分配的内存块指针，失败返回 NULL。
 */
void *debug_calloc(size_t num, size_t size);

/**
 * @brief 调整已分配内存块的大小，并更新记录。
 *
 * @param ptr 要调整的内存块指针。
 * @param size 调整后的大小。
 * @return 调整后的内存块指针，失败返回 NULL。
 */
void *debug_realloc(void *ptr, size_t size);

/**
 * @brief 释放内存块并从记录中移除。
 *
 * @param ptr 要释放的内存块指针。
 */
void debug_free(void *ptr);

/**
 * @brief 检查并打印所有未释放的内存块。
 *
 * 此函数会打印所有泄漏的内存块信息，并计算泄漏的总内存大小。
 * 调用后会清理内存记录链表。
 * @return 泄漏的内存块数量，0 表示无泄漏。
 */
int check_memory_leaks(void);

#ifdef __cplusplus
}
#endif

#endif // DEBUG_MEMORY_H
