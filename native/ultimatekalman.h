#ifndef ULTIMATE_KALMAN
#define ULTIMATE_KALMAN

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifdef BUILD_MEX
#include "mex.h"
#endif

/******************************************************************************/
/* MATRICES                                                                   */
/******************************************************************************/

typedef struct matrix_st {
	int32_t row_dim;
	int32_t col_dim;
	int32_t ld;      // leading dimension
	double* elements;
} matrix_t;

#ifdef NDEBUG
#define matrix_set(A,i,j,v) ( ((A)->elements)[ (j)*((A)->ld) + (i)  ] = (v) )
#define matrix_get(A,i,j)   ( ((A)->elements)[ (j)*((A)->ld) + (i)  ] )
#else
void   matrix_set(matrix_t* A, int32_t i, int32_t j, double v);
double matrix_get(matrix_t* A, int32_t i, int32_t j);
#endif

/*
 * Creates a zero matrix
 */
matrix_t* matrix_create(int32_t rows, int32_t cols);
matrix_t* matrix_create_constant(int32_t rows, int32_t cols, double c);
matrix_t* matrix_create_copy(matrix_t* A);
matrix_t* matrix_create_sub(matrix_t* A, int32_t first_row, int32_t rows, int32_t first_col, int32_t cols);
void      matrix_free(matrix_t* A);

int32_t matrix_rows(matrix_t* A);
int32_t matrix_cols(matrix_t* A);
int32_t matrix_ld  (matrix_t* A);

/*
 * chops the matrix to the given number of rows and columns, from
 * the (0,0) element.
 */
void matrix_mutate_chop(matrix_t* A, int32_t rows, int32_t cols);

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
matrix_t* matrix_create_identity(int32_t rows, int32_t cols);

/*
 * Creates a matrix from a rowwise C matrix.
 */
matrix_t* matrix_create_from_rowwise(double* naked, int32_t rows, int32_t cols);

#ifdef BUILD_MEX
matrix_t* matrix_create_from_mxarray(const mxArray* mx);
mxArray*  matrix_copy_to_mxarray    (matrix_t* A);
#endif


void matrix_print(matrix_t* A, char* format);

matrix_t* matrix_create_vconcat(matrix_t* A, matrix_t* B);

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
static void      farray_free(farray_t* a);
static int64_t   farray_size(farray_t* a);
static int64_t   farray_first_index(farray_t* a);
static int64_t   farray_last_index(farray_t* a);
static void*     farray_get(farray_t* a, int64_t i);
static void*     farray_get_first(farray_t* a);
static void*     farray_get_last(farray_t* a);
static void      farray_append(farray_t* a, void* v);
static void*     farray_drop_first(farray_t* a);
static void*     farray_drop_last(farray_t* a);

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

//typedef struct cov_st {
	//int64_t first; // logical index of first element
//} cov_t;

matrix_t* cov_weigh(matrix_t* cov, char cov_type, matrix_t* A);
//matrix_t* cov_weigh(cov_t* cov, matrix_t* A);

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

typedef struct step_st {
	int64_t step; // logical step number
	int32_t dimension;

	matrix_t* Rdiag;
	matrix_t* Rsupdiag;
	matrix_t* y;

	matrix_t* Rbar;
	matrix_t* ybar;

	matrix_t* state;
	matrix_t* covariance;
} step_t;

step_t* step_create();
void    step_free(step_t* step);

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

typedef struct kalman_st {
	//int64_t first; // logical index of first element
	farray_t* steps;
	step_t*   current;
} kalman_t;

kalman_t* kalman_create    ();
void      kalman_free      (kalman_t* kalman);

int64_t   kalman_earliest  (kalman_t* kalman);
int64_t   kalman_latest    (kalman_t* kalman);

void      kalman_evolve    (kalman_t* kalman, int32_t n_i, matrix_t* H_i, matrix_t* F_i, matrix_t* c_i, matrix_t* K_i, char K_type);
void      kalman_observe   (kalman_t* kalman, matrix_t* G_i, matrix_t* o_i, matrix_t* C_i, char C_type);
void      kalman_smooth    (kalman_t* kalman);
matrix_t* kalman_estimate  (kalman_t* kalman, int64_t si);
matrix_t* kalman_covariance(kalman_t* kalman, int64_t si);
void      kalman_forget    (kalman_t* kalman, int64_t si);
void      kalman_rollback  (kalman_t* kalman, int64_t si);
matrix_t* kalman_perftest(kalman_t* kalman,
		                      matrix_t* H, matrix_t* F, matrix_t* c, matrix_t* K, char K_type,
		                      matrix_t* G, matrix_t* o,              matrix_t* C, char C_type,
									        int32_t count, int32_t decimation);


/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
