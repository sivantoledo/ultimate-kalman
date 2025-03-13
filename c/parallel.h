#ifndef PARALLEL_H
#define PARALLEL_H

#include "concurrent_set.h"

void parallel_set_thread_limit(int number_of_threads);
void parallel_set_blocksize   (int blocksize_in);

void foreach_in_range    (void (*func)(void*,        int, size_t, size_t), void* array ,               int length, size_t n);
void foreach_in_range_two(void (*func)(void*, void*, int, size_t, size_t), void* array1, void* array2, int length, size_t n);

void prefix_sums_pointers(void* (*f)(void*, void*), void** input, void** sums, concurrent_set_t* create_array , int length, int stride);

// Opaque pointer for clients
struct spin_mutex_st;
typedef struct spin_mutex_st spin_mutex_t;

spin_mutex_t* spin_mutex_create ();
void          spin_mutex_lock   (spin_mutex_t* mutex);
void          spin_mutex_unlock (spin_mutex_t* mutex);
void          spin_mutex_destroy(spin_mutex_t* mutex);

#endif
