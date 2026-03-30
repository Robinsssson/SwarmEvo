#include "alg_memalloc.h"
#include <stdlib.h>
static alg_memalloc_hook g_memalloc_handle = {.alg_memalloc_calloc = calloc,
                                              .alg_memalloc_free = free,
                                              .alg_memalloc_malloc = malloc,
                                              .alg_memalloc_realloc = realloc};

void alg_memalloc_init(alg_memalloc_hook *handle) {
    if (handle == NULL) {
        g_memalloc_handle.alg_memalloc_calloc = calloc;
        g_memalloc_handle.alg_memalloc_free = free;
        g_memalloc_handle.alg_memalloc_malloc = malloc;
        g_memalloc_handle.alg_memalloc_realloc = realloc;
    } else {
        g_memalloc_handle.alg_memalloc_realloc = handle->alg_memalloc_realloc;
        g_memalloc_handle.alg_memalloc_calloc = handle->alg_memalloc_calloc;
        g_memalloc_handle.alg_memalloc_free = handle->alg_memalloc_free;
        g_memalloc_handle.alg_memalloc_malloc = handle->alg_memalloc_malloc;
    }
}

void *alg_malloc(size_t size) {
    return g_memalloc_handle.alg_memalloc_malloc(size);
}
void alg_free(void *ptr) {
    g_memalloc_handle.alg_memalloc_free(ptr);
}
void *alg_realloc(void *ptr, size_t size) {
    return g_memalloc_handle.alg_memalloc_realloc(ptr, size);
}
void *alg_calloc(size_t number, size_t size) {
    return g_memalloc_handle.alg_memalloc_calloc(number, size);
}
