#ifndef FLEXIBLE_ARRAYS_H
#define FLEXIBLE_ARRAYS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef BUILD_MEX
#include "mex.h"
#endif

/******************************************************************************/
/* FLEXIBLE ARRAYS OF POINTERS                                                */
/******************************************************************************/

typedef struct farray_st {
	int64_t first; // logical index of first element
	int64_t start; // physical index
	int64_t end;   // physical index
	int64_t array_size;
	void**  elements;
} farray_t;

farray_t* farray_create();
void      farray_free(farray_t* a);
int64_t   farray_size(farray_t* a);
int64_t   farray_first_index(farray_t* a);
int64_t   farray_last_index(farray_t* a);
void*     farray_get(farray_t* a, int64_t i);
void*     farray_get_first(farray_t* a);
void*     farray_get_last(farray_t* a);
void      farray_append(farray_t* a, void* v);
void*     farray_drop_first(farray_t* a);
void*     farray_drop_last(farray_t* a);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
