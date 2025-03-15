/*
 * flexible_arrays.c
 *
 * (C) Sivan Toledo, 2022-2025
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
// for "unused" attribute
#define __attribute__(x)
#include <float.h>
// string.h for memcpy
#else
#include <unistd.h>
#endif

#ifdef BUILD_MEX
#include "mex.h"

static char assert_msg[128];
static void mex_assert(int c, int line) {
    if (!c) {
        sprintf(assert_msg,"Assert failed in %s line %d",__FILE__,line);
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:assertion",assert_msg);
    }
}

#define assert(c) mex_assert((c),__LINE__)
#else
#include <assert.h>
#endif

#include "flexible_arrays.h"

/******************************************************************************/
/* FLEXIBLE ARRAYS OF POINTERS                                                */
/******************************************************************************/

#define FARRAY_INITIAL_SIZE 1024

farray_t* farray_create() {
  farray_t *a = malloc(sizeof(farray_t));
  assert(a != NULL);
  a->array_size = FARRAY_INITIAL_SIZE;
  a->first = 0; // once we append, the first element will be 0
  a->start = 0; // in preparation for elements
  a->end = -1;
  a->elements = calloc(a->array_size, sizeof(void*));
  assert(a->elements != NULL);
  return a;
}

void farray_free(farray_t *a) {
  free(a->elements);
  free(a);
}

farray_index_t farray_size(farray_t *a) {
  if (a->end < a->start)
    return 0;
  return (a->end) - (a->start) + 1;
}

farray_index_t farray_first_index(farray_t *a) {
  return (a->first);
}

farray_index_t farray_last_index(farray_t *a) {
  return (a->first) + farray_size(a) - 1;
}

void* farray_get(farray_t *a, farray_index_t i) {

  assert(i >= a->first);
  assert(i < (a->first) + farray_size(a));

  farray_index_t offset = i - (a->first);
  farray_index_t physical = (a->start) + offset;

  return (a->elements)[physical];
}

void* farray_get_first(farray_t *a) {
  return farray_get(a, a->first);
}

void* farray_get_last(farray_t *a) {
  return farray_get(a, (a->first) + farray_size(a) - 1);
}

void farray_append(farray_t *a, void *v) {
  farray_index_t i;

  if (a->end >= (a->array_size) - 1) {
    // shift back or extend

    farray_index_t logical_size = farray_size(a);
    farray_index_t array_size = a->array_size;

    // if (debug) printf("append l=%lld p=%lld\n",logical_size,array_size);

    if (logical_size <= array_size / 2) { // can just shift back less than half the array
      for (i = 0; i < logical_size; i++) {
        (a->elements)[i] = (a->elements)[(a->start) + i];
      }
      a->start = 0;
      a->end = logical_size - 1;
      //if (debug) printf("farray shifted back, physical size %lld pointers, logical size %lld\n",a->array_size,logical_size);
    } else { // array is currently more than half full, realloc at a larger size
      a->array_size *= 2;
      a->elements = realloc(a->elements, (a->array_size) * sizeof(void*));
      assert(a->elements != NULL);
      //if (debug) printf("farray reallocated, new physical size is %lld pointers, logical size %lld\n",a->array_size,logical_size);
    }
  }

  assert(a->end < (a->array_size) - 1);

  (a->end)++;
  (a->elements)[a->end] = v;
}

void* farray_drop_first(farray_t *a) {
  void *r = (a->elements)[a->start];
  (a->first)++;
  (a->start)++;
  return r;
}

void* farray_drop_last(farray_t *a) {
  void *r = (a->elements)[a->end];
  (a->end)--;
  return r;
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
