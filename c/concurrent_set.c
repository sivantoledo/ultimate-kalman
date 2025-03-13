/*
 * A simple concurrent set data structure, to keep track of objects
 * that are created during the parallel prefix sum operation (the elements
 * are structures, not value types) so that we can release them at the end of
 * the operation.
 *
 * The current code may fail if the representation uses an array whose indexes
 * cannot be represented by 32-bit integers, because it will try to put all
 * the elements in the first 2 billion or so cells.
 *
 * Copyright (c) Sivan Toledo and Shahaf Gargir 2024-2025
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "parallel.h"

typedef struct concurrent_set_st {
  parallel_index_t size;
  void**           pointers;
  spin_mutex_t**   locks; // pointers to locks, to allow locks of any type
  void             (*foreach)(void*);
} concurrent_set_t;

/*
 * FNV-1a hash of an address, to generate a random integer
 * 
 * Based on code suggested by Grok 3. 
 */

static uint32_t hash_uint32(uint32_t value) {
  const uint32_t FNV_PRIME = 16777619u;
  const uint32_t FNV_OFFSET = 3141592653u; // Different offset for distinctness

  uint32_t hash = FNV_OFFSET;
  for (size_t i = 0; i < sizeof(uint32_t); i++) {
    hash ^= (value & 0xFF);
    hash *= FNV_PRIME;
    value >>= 8;
  }
  return hash;
}

static void concurrent_set_parallel_init(void *set_v, parallel_index_t length, parallel_index_t start, parallel_index_t end) {
  concurrent_set_t *set = (concurrent_set_t*) set_v;
  for (parallel_index_t i = start; i < end; i++) {
    (set->pointers)[i] = NULL;
    (set->locks)[i] = spin_mutex_create();
  }
}

static void concurrent_set_parallel_destroy(void *set_v, parallel_index_t length, parallel_index_t start, parallel_index_t end) {
  concurrent_set_t *set = (concurrent_set_t*) set_v;
  for (parallel_index_t i = start; i < end; i++) {
    spin_mutex_destroy((set->locks)[i]);
  }
}

static void concurrent_set_parallel_foreach(void *set_v, parallel_index_t length, parallel_index_t start, parallel_index_t end) {
  concurrent_set_t *set = (concurrent_set_t*) set_v;
  for (parallel_index_t i = start; i < end; i++) {
    if ((set->pointers)[i] != NULL) {
      (*(set->foreach))((set->pointers)[i]);
    }
  }
}

concurrent_set_t* concurrent_set_create(parallel_index_t capacity, void (*foreach)(void*)) {
  concurrent_set_t *set = (concurrent_set_t*) malloc(sizeof(concurrent_set_t));
  set->size = capacity * 10; // expansion to reduce contention
  set->foreach = foreach;
  set->pointers = (void**) malloc((set->size) * sizeof(void*));
  set->locks = (spin_mutex_t**) malloc((set->size) * sizeof(spin_mutex_t*));

  //parallel_for_c(la, NULL, 0, k, BLOCKSIZE, parallelInit);
  foreach_in_range(concurrent_set_parallel_init, set, set->size, set->size);

  return set;
}

void concurrent_set_free(concurrent_set_t *set) {
  foreach_in_range(concurrent_set_parallel_destroy, set, set->size, set->size);
  //parallel_for_c(la, NULL, 0, la->rows, BLOCKSIZE, parallelDestroy);
  free(set->locks);
  free(set->pointers);
  free(set);
}

void concurrent_set_insert(concurrent_set_t *set, void *element) {
  uint32_t inserted = 0;
  uint32_t h = (uint32_t) (uintptr_t) element;
  uint32_t i;
  do {
    h = hash_uint32(h);
    i = h % (set->size);

    spin_mutex_lock((set->locks)[i]);
    if ((set->pointers)[i] == NULL) {
      (set->pointers)[i] = element;
      inserted = 1;
    }
    spin_mutex_unlock((set->locks)[i]);
  } while (!inserted);
}

void concurrent_set_foreach(concurrent_set_t *set) {
  foreach_in_range(concurrent_set_parallel_foreach, set, set->size, set->size);
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
