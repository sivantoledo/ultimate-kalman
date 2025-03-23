/*
 * parallel.h
 *
 * Definitions of the parallel primitives for a collection of Kalman
 * filters and smoothers.
 *
 * Copyright (c) 2024-2025 Sivan Toledo and Shahaf Gargir
 */

#ifndef PARALLEL_H
#define PARALLEL_H

#ifdef PARALLEL_TYPE_INT32
typedef int32_t parallel_index_t;
#endif

#ifdef PARALLEL_INDEX_TYPE_UINT32
typedef uint32_t parallel_index_t;
#endif

#ifdef PARALLEL_INDEX_TYPE_INT64
typedef int64_t parallel_index_t;
#endif

#ifdef PARALLEL_INDEX_TYPE_UINT64
typedef uint64_t parallel_index_t;
#endif

#include "concurrent_set.h"

void parallel_set_thread_limit(int number_of_threads);
void parallel_set_blocksize   (int blocksize_in);

void foreach_in_range    (void (*func)(void*,        parallel_index_t, parallel_index_t, parallel_index_t), void* array ,               parallel_index_t length, parallel_index_t n);
void foreach_in_range_two(void (*func)(void*, void*, parallel_index_t, parallel_index_t, parallel_index_t), void* array1, void* array2, parallel_index_t length, parallel_index_t n);

void prefix_sums_pointers(void* (*f)(void*, void*), void** input, void** sums, concurrent_set_t* create_array , parallel_index_t length, int stride);

// Opaque pointer for clients
struct spin_mutex_st;
typedef struct spin_mutex_st spin_mutex_t;

spin_mutex_t* spin_mutex_create ();
void          spin_mutex_lock   (spin_mutex_t* mutex);
void          spin_mutex_unlock (spin_mutex_t* mutex);
void          spin_mutex_destroy(spin_mutex_t* mutex);

#endif
