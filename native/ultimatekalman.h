#ifndef ULTIMATE_KALMAN
#define ULTIMATE_KALMAN

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

/******************************************************************************/
/* MATRICES                                                                   */
/******************************************************************************/

typedef struct matrix_st {
	int32_t row_dim;
	int32_t col_dim;
	double* elements;
} matrix_t;

#ifdef NDEBUG
#define matrix_set(A,i,j,v) ( ((A)->elements)[ (j)*((A)->row_dim) + (i)  ] = (v) )
#define matrix_get(A,i,j)   ( ((A)->elements)[ (j)*((A)->row_dim) + (i)  ] )
#else
void   matrix_set(matrix_t* A, int32_t i, int32_t j, double v);
double matrix_get(matrix_t* A, int32_t i, int32_t j);
#endif

/*
 * Creates a zero matrix
 */
matrix_t* matrix_create(int32_t rows, int32_t cols);

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
matrix_t* matrix_create_identity(int32_t rows, int32_t cols);

/*
 * Creates a matrix from a rowwise C matrix.
 */
matrix_t* matrix_create_from_rowwise(double* naked, int32_t rows, int32_t cols);

void matrix_print(matrix_t* A, FILE* f, char* format);

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

static farray_t* farray_create();
static int64_t   farray_size(farray_t* a);
static int64_t   farray_first_index(farray_t* a);
static int64_t   farray_last_index(farray_t* a);
static void*     farray_get(farray_t* a, int64_t i);
static void*     farray_get_first(farray_t* a);
static void*     farray_get_last(farray_t* a);
static void      farray_append(farray_t* a, void* v);
static void      farray_drop_first(farray_t* a);

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

typedef struct cov_st {
	int64_t first; // logical index of first element
	//farray_t* steps;
} cov_t;

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

typedef struct step_st {
	int64_t first; // logical index of first element
	//farray_t* steps;
} step_t;

static step_t* step_create();

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

typedef struct kalman_st {
	//int64_t first; // logical index of first element
	farray_t* steps;
	step_t*   current;
} kalman_t;

kalman_t* kalman_create();
void      advance(kalman_t* kalman, int32_t dim, double timestamp);
void      evolve(kalman_t* kalman, matrix_t* H, matrix_t* F, matrix_t* be, cov_t* Ce);
void      evolve_simple(kalman_t* kalman, matrix_t* F, matrix_t* be, cov_t* Ce);
void      observe_nothing(kalman_t* kalman);
void      observe(kalman_t* kalman, matrix_t* G, matrix_t* bo, cov_t* Co);
matrix_t* filter(kalman_t* kalman);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
