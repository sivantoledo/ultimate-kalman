
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <blas.h>
#include <lapack.h>

#include "ultimatekalman.h"

/******************************************************************************/
/* UTILITIES                                                                  */
/******************************************************************************/

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

/******************************************************************************/
/* MATRICES                                                                   */
/******************************************************************************/

#ifdef NDEBUG
// defined in the header file ...
#else
static void matrix_set(matrix_t* A, int32_t i, int32_t j, double v) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < A->row_dim);
	assert(j < A->col_dim);
	(A->elements)[ j*(A->row_dim) + i ] = v;
}

static double matrix_get(matrix_t* A, int32_t i, int32_t j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < A->row_dim);
	assert(j < A->col_dim);
	return (A->elements)[ j*(A->row_dim) + i ];
}
#endif

/*
 * Creates a zero matrix
 */
static matrix_t* matrix_create(int32_t rows, int32_t cols) {
	matrix_t* A = malloc(sizeof(matrix_t));
	assert( A!= NULL );
	A->row_dim = rows;
	A->col_dim = cols;
	A->elements = calloc(rows*cols, sizeof(double));
	assert( A->elements != NULL );
	return A;
}

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
static matrix_t* matrix_create_identity(int32_t rows, int32_t cols) {
	int32_t i,j;

	matrix_t* I = matrix_create(rows,cols);

	for (i=0; i<MIN(rows,cols); i++) {
		matrix_set(I,i,i,1.0);
	}

	return I;
}

static matrix_t* matrix_create_from_rowwise(double* naked, int32_t rows, int32_t cols) {
	int32_t i,j;

	matrix_t* A = malloc(sizeof(matrix_t));
	assert( A!= NULL );
	A->row_dim = rows;
	A->col_dim = cols;
	A->elements = calloc(rows*cols, sizeof(double));
	assert( A->elements != NULL );

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(A,i,j, naked[ i*cols + j ]);
		}
	}

	return A;
}


static void matrix_print(matrix_t* A, FILE* f, char* format) {
	int32_t i,j;

	if (f == NULL)      f      = stdout;
	if (format == NULL) format = "%f";

	for (i=0; i<A->row_dim; i++) {
		for (j=0; j<A->col_dim; j++) {
			fprintf(f,format,matrix_get(A,i,j));
			fprintf(f," ");
		}
		fprintf(f,"\n");
	}
}

/******************************************************************************/
/* FLEXIBLE ARRAYS OF POINTERS                                                */
/******************************************************************************/

#define FARRAY_INITIAL_SIZE 1024

static farray_t* farray_create() {
	farray_t* a = malloc(sizeof(farray_t));
	assert( a != NULL );
	a->array_size = FARRAY_INITIAL_SIZE;
	a->first      =  0; // once we append, the first element will be 0
	a->start      =  0; // in preparation for elements
	a->end        = -1;
	a->elements   = calloc(a->array_size, sizeof(void*));
	assert( a->elements != NULL );
	return a;
}

static int64_t farray_size(farray_t* a) {
	if (a->end < a->start) return 0;
	return (a->end) - (a->start) + 1;
}

static int64_t farray_first_index(farray_t* a) {
	return (a->first);
}

static int64_t farray_last_index(farray_t* a) {
	return (a->first) + farray_size(a) -1;
}

static void* farray_get(farray_t* a, int64_t i) {

	assert( i >= a->first );
	assert( i < (a->first) + farray_size(a) );

	int64_t offset = i - (a->first);
	int64_t physical = (a->start) + offset;

	return (a->elements)[ physical ];
}

static void* farray_get_first(farray_t* a) {
	return farray_get(a, a->first);
}

static void* farray_get_last(farray_t* a) {
	return farray_get(a, (a->first) + farray_size(a) -1 );
}

static void farray_append(farray_t* a, void* v) {
	int64_t i;

	if (a->end >= (a->array_size)-1) {
		// shift back or extend

		int64_t logical_size = farray_size(a);
		int64_t array_size   = a->array_size;

		// printf("append l=%lld p=%lld\n",logical_size,array_size);

		if (logical_size <= array_size/2) { // can just shift back less than half the array
			for (i=0; i<logical_size; i++) {
				(a->elements)[ i ] = (a->elements)[ (a->start) + i ];
			}
			a->start = 0;
			a->end   = logical_size-1;
			//printf("farray shifted back, physical size %lld pointers, logical size %lld\n",a->array_size,logical_size);
		} else { // array is currently more than half full, realloc at a larger size
			a->array_size *= 2;
			a->elements = realloc( a->elements, (a->array_size)*sizeof(void*) );
			assert( a->elements != NULL);
			//printf("farray reallocated, new physical size is %lld pointers, logical size %lld\n",a->array_size,logical_size);
		}
	}

	assert( a->end < (a->array_size)-1 );

	(a->end)++;
	(a->elements)[ a->end ] = v;
}

static void farray_drop_first(farray_t* a) {
	(a->first)++;
	(a->start)++;
}

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

static step_t* step_create() {
	step_t* s = malloc(sizeof(step_t));
	assert( s != NULL );
	return s;
}

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

kalman_t* kalman_create() {
	kalman_t* kalman = malloc(sizeof(kalman_t));
	assert( kalman != NULL );
	kalman->steps   = farray_create();
	kalman->current = NULL;
	return kalman;
}

void kalman_free(kalman_t* kalman) {
	printf("waning: kalman_free not yet implemented\n");
	//if (kalman->steps != NULL)
}

void advance(kalman_t* kalman, int32_t dim, double timestamp) {
	kalman->current = step_create(dim, timestamp);
}

void evolve(kalman_t* kalman, matrix_t* H, matrix_t* F, matrix_t* be, cov_t* Ce) {
}

void evolve_simple(kalman_t* kalman, matrix_t* F, matrix_t* be, cov_t* Ce) {
}

void observe_nothing(kalman_t* kalman) {
}

void observe(kalman_t* kalman, matrix_t* G, matrix_t* bo, cov_t* Co) {
}

matrix_t* filter(kalman_t* kalman) {
}

/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

int main(int argc, char *argv[]) {
  printf("hello world\n");

#if 0
  matrix_t* A = matrix_create(7,8);
  matrix_t* I = matrix_create_identity(7,8);
  matrix_set(I,6,7,13.0);
  matrix_print(I,NULL,NULL);
  //matrix_set(I,7,0,13.0);
#endif

  farray_t* a = farray_create();
  printf("farray size %lld\n",farray_size(a));

  farray_append(a,(void*) 0);
  farray_append(a,(void*) 1);
  farray_append(a,(void*) 2);
  farray_append(a,(void*) 3);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));
  for (int i=farray_first_index(a); i<=farray_last_index(a); i++) {
  	printf("farray get(%d) %lld\n",i,(int64_t) farray_get(a,i));
  }

  farray_drop_first(a);
  farray_drop_first(a);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));
  for (int i=farray_first_index(a); i<=farray_last_index(a); i++) {
  	printf("farray get(%d) %lld\n",i,(int64_t) farray_get(a,i));
  }

	printf("farray get_first %lld\n",(int64_t) farray_get_first(a));
	printf("farray get_last  %lld\n",(int64_t) farray_get_last(a));

  printf("  farray append 10\n");
  for (int i=0; i<10; i++) farray_append(a,(void*)77);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));

  printf("  farray drop 9\n");
  for (int i=0; i<9; i++) farray_drop_first(a);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));

  printf("  farray append 128\n");
  for (int i=0; i<128; i++) farray_append(a,(void*)77);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));

  double A[] = {
    1, 2, 3,
  	4, 5, 6
  };

  matrix_t* W = matrix_create_from_rowwise(A,2,3);
  matrix_print(W,NULL,NULL);
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

