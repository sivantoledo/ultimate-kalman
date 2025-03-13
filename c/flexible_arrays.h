#ifndef FLEXIBLE_ARRAYS_H
#define FLEXIBLE_ARRAYS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/******************************************************************************/
/* FLEXIBLE ARRAYS OF POINTERS                                                */
/******************************************************************************/

typedef struct farray_st {
	farray_index_t first; // logical index of first element
	farray_index_t start; // physical index
	farray_index_t end;   // physical index
	farray_index_t array_size;
	void**  elements;
} farray_t;

farray_t* farray_create();
void      farray_free(farray_t* a);
farray_index_t   farray_size(farray_t* a);
farray_index_t   farray_first_index(farray_t* a);
farray_index_t   farray_last_index(farray_t* a);
void*     farray_get(farray_t* a, farray_index_t i);
void*     farray_get_first(farray_t* a);
void*     farray_get_last(farray_t* a);
void      farray_append(farray_t* a, void* v);
void*     farray_drop_first(farray_t* a);
void*     farray_drop_last(farray_t* a);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
