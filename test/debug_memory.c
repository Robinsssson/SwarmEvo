#include "debug_memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// 定义节点结构以跟踪已分配的内存块
typedef struct MemNode {
    void *ptr;            // 内存块指针
    size_t size;          // 分配大小
    struct MemNode *next; // 下一个节点
} MemNode;

// 全局链表头节点
static MemNode *head = NULL;

// 添加记录到链表
static void add_mem_record(void *ptr, size_t size) {
    MemNode *node = (MemNode *)malloc(sizeof(MemNode));
    node->ptr = ptr;
    node->size = size;
    node->next = head;
    head = node;
}

// 从链表中移除记录
static void remove_mem_record(void *ptr) {
    MemNode *prev = NULL, *current = head;

    while (current) {
        if (current->ptr == ptr) {
            if (prev) {
                prev->next = current->next;
            } else {
                head = current->next;
            }
            free(current);
            return;
        }
        prev = current;
        current = current->next;
    }
}

// 自定义 malloc 函数
void *debug_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr) {
        add_mem_record(ptr, size);
    }
    return ptr;
}

// 自定义 calloc 函数
void *debug_calloc(size_t num, size_t size) {
    void *ptr = calloc(num, size);
    if (ptr) {
        add_mem_record(ptr, num * size);
    }
    return ptr;
}

// 自定义 realloc 函数
void *debug_realloc(void *ptr, size_t size) {
    if (ptr) {
        remove_mem_record(ptr);
    }
    void *new_ptr = realloc(ptr, size);
    if (new_ptr) {
        add_mem_record(new_ptr, size);
    }
    return new_ptr;
}

// 自定义 free 函数
void debug_free(void *ptr) {
    if (ptr) {
        remove_mem_record(ptr);
        free(ptr);
    }
}

// 清理所有记录并检查是否有未释放的内存，返回泄漏块数
int check_memory_leaks(void) {
    MemNode *current = head;
    int leaked = 0;

    if (!current) {
        printf("No memory leaks detected.\n");
        return 0;
    }

    printf("Memory leaks detected:\n");
    size_t total_leaked = 0;
    while (current) {
        printf("Leaked block at %p, size: %zu bytes\n", current->ptr, current->size);
        total_leaked += current->size;
        leaked++;
        current = current->next;
    }
    printf("Total leaked memory: %zu bytes\n", total_leaked);

    // 清理链表
    current = head;
    while (current) {
        MemNode *to_free = current;
        current = current->next;
        free(to_free);
    }
    head = NULL;
    return leaked;
}
