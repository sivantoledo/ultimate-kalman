#ifndef CONCURRENT_SET_H
#define CONCURRENT_SET_H

#include <stdint.h>
#include <stdlib.h>
 
/*
 * These definitions of the integer data type are from parallel.h, but
 * we can't include it because it would be circular.
 */

#ifdef PARALLEL_INDEX_TYPE_INT32
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

struct concurrent_set_st;
typedef struct concurrent_set_st concurrent_set_t;

concurrent_set_t* concurrent_set_create (parallel_index_t capacity, void (*foreach)(void*));
void              concurrent_set_free   (concurrent_set_t* set);
void              concurrent_set_insert (concurrent_set_t* set, void* element);
void              concurrent_set_foreach(concurrent_set_t* set);

#endif
