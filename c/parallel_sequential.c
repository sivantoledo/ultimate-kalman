/*
 * parallel_sequential.c
 *
 * A sequential implementation of the parallel primitives (mostly for computing speedups).
 *
 * Copyright (c) 2024-2025 Sivan Toledo and Shahaf Gargir
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "parallel.h"

void parallel_set_thread_limit(int number_of_threads) {
}
void parallel_set_blocksize(int blocksize_in) {
}

void foreach_in_range(void (*f)(void*, parallel_index_t, parallel_index_t, parallel_index_t),
                      void *array, parallel_index_t length, parallel_index_t n) {
  (*f)(array, length, 0, n);
}

void foreach_in_range_two(void (*f)(void*, void*, parallel_index_t, parallel_index_t, parallel_index_t),
                          void *array1, void *array2, parallel_index_t length, parallel_index_t n) {
  (*f)(array1, array2, length, 0, n);
}

void prefix_sums_pointers(void* (*f)(void*, void*),
                          void **input, void **sums,
                          concurrent_set_t *created_elements,
                          parallel_index_t length, int stride) {
  parallel_index_t i, j;
  void *sum = NULL; // neutral element when operating on pointers

  for (i = 0; i < length; i++) {
    if (stride == 1) {
      j = i;
    } else {
      j = length - 1 - i;
    }
    //fprintf(stderr,">>> %d %d\n",i,j);
    void *temp = f(sum, input[j]);
    if ((sum != NULL) && (input[j] != NULL))
      concurrent_set_insert(created_elements, temp); // the first element is combined with NULL so f returns it, not a new element
    sums[i] = temp;
    sum = temp;
  }
}

spin_mutex_t* spin_mutex_create()            { return NULL; }
void spin_mutex_lock(spin_mutex_t *mutex)    {}
void spin_mutex_unlock(spin_mutex_t *mutex)  {}
void spin_mutex_destroy(spin_mutex_t *mutex) {}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
