
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <blas.h>
#include <lapack.h>

#include "ultimatekalman.h"

#ifdef BUILD_MEX
#include "mex.h"

static char assert_msg[81];
static void mex_assert(int c) {
	if (!c) {
		sprintf(assert_msg,"Assert failed in %s line %d",__FILE__,__LINE__);
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:assertion",assert_msg);
	}
}

#define assert(c) mex_assert((c))
#define blas_int_t mwSignedIndex
#else
#define blas_int_t int
#endif


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
void matrix_set(matrix_t* A, int32_t i, int32_t j, double v) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < A->row_dim);
	assert(j < A->col_dim);
	(A->elements)[ j*(A->row_dim) + i ] = v;
}

double matrix_get(matrix_t* A, int32_t i, int32_t j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < A->row_dim);
	assert(j < A->col_dim);
	return (A->elements)[ j*(A->row_dim) + i ];
}
#endif

/*
 * Creates a matrix with undefined elements
 */
matrix_t* matrix_create(int32_t rows, int32_t cols) {
	matrix_t* A = malloc(sizeof(matrix_t));
	assert( A!= NULL );
	A->row_dim = rows;
	A->col_dim = cols;
	A->elements = malloc(rows*cols*sizeof(double));
	assert( A->elements != NULL );
	return A;
}

void matrix_free(matrix_t* A) {
	if (A==NULL) return;
	free( A->elements );
	free( A );
}

int32_t matrix_rows(matrix_t* A) { return A->row_dim; }
int32_t matrix_cols(matrix_t* A) { return A->col_dim; }

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
matrix_t* matrix_create_identity(int32_t rows, int32_t cols) {
	int32_t i;

	matrix_t* I = matrix_create(rows,cols);

	for (i=0; i<MIN(rows,cols); i++) {
		matrix_set(I,i,i,1.0);
	}

	return I;
}

matrix_t* matrix_create_from_rowwise(double* naked, int32_t rows, int32_t cols) {
	int32_t i,j;

	matrix_t* A = matrix_create(rows,cols);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(A,i,j, naked[ i*cols + j ]);
		}
	}

	return A;
}

matrix_t* matrix_create_constant(int32_t rows, int32_t cols, double c) {
	int32_t i,j;

	matrix_t* A = matrix_create(rows,cols);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(A,i,j, c);
		}
	}

	return A;
}

void matrix_print(matrix_t* A, char* format) {
	int32_t i,j;

	printf("matrix_print %d %d\n",A->row_dim,A->col_dim);

	if (format == NULL) format = "%f";

	for (i=0; i<A->row_dim; i++) {
		for (j=0; j<A->col_dim; j++) {
			printf(format,matrix_get(A,i,j));
			printf(" ");
		}
		printf("\n");
	}
}

matrix_t* matrix_create_copy(matrix_t* A) {
	int i,j;

	int rows  = (A->row_dim);
	int cols  = (A->col_dim);

	matrix_t* C = matrix_create(rows,cols);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(A,i,j));
		}
	}

	return C;
}

matrix_t* matrix_create_vconcat(matrix_t* A, matrix_t* B) {
	int i,j;

	assert( A->col_dim == B->col_dim );

	int rows  = (A->row_dim) + (B->row_dim);
	int cols  = (A->col_dim);
	int Arows = A->row_dim;

	matrix_t* C = matrix_create(rows,cols);

	for (i=0; i<Arows; i++) {
		for (j=0; j<C->col_dim; j++) {
			matrix_set(C,i,j,matrix_get(A,i,j));
		}
	}

	for (   ; i<C->row_dim; i++) {
		for (j=0; j<C->col_dim; j++) {
			matrix_set(C,i,j,matrix_get(B,i-Arows,j));
		}
	}

	return C;
}


#ifdef BUILD_MEX
matrix_t* matrix_create_from_mxarray(const mxArray* mx) {
	int32_t rows = (int32_t) mxGetM(mx);
	int32_t cols = (int32_t) mxGetN(mx);

	if (rows==0 || cols==0) return NULL;

	printf("== create from mx %d %d\n",rows,cols);

	matrix_t* A = malloc(sizeof(matrix_t));
	assert( A!= NULL );
	A->row_dim = rows;
	A->col_dim = cols;
	A->elements = calloc(rows*cols, sizeof(double));
	assert( A->elements != NULL );

	memcpy( A->elements, mxGetPr(mx), rows*cols*sizeof(double));

	return A;
}
#endif


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

static void farray_free(farray_t* a) {
  free( a->elements );
  free( a );
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

matrix_t* cov_weigh(cov_t* cov, matrix_t* A) {
	printf("cov_weigh not implemented yet\n");
	assert(false);
	return NULL;
}

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
	printf("waning: kalman_free not fully implemented yet (steps not processed)\n");
	farray_free( kalman->steps );
	// step_free( kalman->current );
	free( kalman );
}

int64_t   kalman_earliest(kalman_t* kalman) {
	if ( farray_size(kalman->steps) == 0 ) return -1;
	step_t* s = farray_get_first(kalman->steps);
	return s->step;
}

int64_t   kalman_latest(kalman_t* kalman) {
	if ( farray_size(kalman->steps) == 0 ) return -1;
	step_t* s = farray_get_last(kalman->steps);
	return s->step;
}

void kalman_evolve(kalman_t* kalman, int32_t n_i, matrix_t* H_i, matrix_t* F_i, matrix_t* c_i, cov_t* K_i) {
	printf("kalman_evolve n_i = %d\n",n_i);
	printf("kalman_evolve H_i = %08x\n",H_i);
	printf("kalman_evolve F_i = %08x\n",F_i);
	printf("kalman_evolve i_i = %08x\n",c_i);
	printf("kalman_evolve K_i = %08x\n",K_i);

	kalman->current = step_create();

	if (farray_size(kalman->steps)==0) {
		printf("kalman_evolve first step\n");
		kalman->current->step = 0;
		return;
	}

	step_t* imo = farray_get_last( kalman->steps );
	kalman->current->step = (imo->step) + 1;

	printf("kalman_evolve step = %d\n",kalman->current->step);

	assert(H_i!=NULL);
	assert(F_i!=NULL);
	assert(c_i!=NULL);
	assert(K_i!=NULL);

	printf("kalman_evolve F_i = %08x %d %d\n",F_i,F_i->row_dim,F_i->col_dim);
	matrix_print(F_i,NULL);
return;
	matrix_t* V_i_H_i = cov_weigh(K_i,H_i);
	matrix_t* V_i_F_i = cov_weigh(K_i,F_i);
	matrix_t* V_i_c_i = cov_weigh(K_i,c_i);

	matrix_t* A;
	matrix_t* B;
	matrix_t* y;

	if ( imo->Rdiag != NULL ) {
		int32_t z_i = matrix_rows( imo->Rdiag );
		A = matrix_create_vconcat( imo->Rdiag                             , V_i_F_i );
		B = matrix_create_vconcat( matrix_create_constant( z_i, n_i, 0.0 ), V_i_H_i );
		y = matrix_create_vconcat( imo->y                                 , V_i_c_i );
	} else {
		A = matrix_create_copy( V_i_F_i );
		B = matrix_create_copy( V_i_H_i );
		y = matrix_create_copy( V_i_c_i );
	}

	/*
	 * QR factorization
	 */

	blas_int_t M,N,K,LDA,LDC,LWORK,INFO;

	M   = matrix_rows(A);
	N   = matrix_cols(A);
	LDA = matrix_rows(A);

	matrix_t* TAU = matrix_create(M,1);

	double WORK_SCALAR;
	LWORK = -1; // tell lapack to compute the size of the work area required

#if 0
	dgeqrf(&M, &N, A->elements, &LDA, TAU->elements, &WORK_SCALAR, &LWORK, &INFO);

	if (INFO != 0) printf("dgeqrf INFO=%d\n",INFO);
	assert(INFO==0);

	printf("dgeqrf requires %f words in WORK\n",WORK_SCALAR);

	LWORK = (int32_t) WORK_SCALAR;
	matrix_t* WORK = matrix_create(LWORK,1);

	dgeqrf(&M, &N, A->elements, &LDA, TAU->elements, WORK->elements, &LWORK, &INFO);

	matrix_free(WORK);

	if (INFO != 0) printf("dgeqrf INFO=%d\n",INFO);
	assert(INFO==0);
#endif
	// left (pre) multiplication by Q^T
	//dormqr("L", "T", &M, &N, &K, A->elements, &LDA, tau, B->elements, &LDC, &WORK, &LWORD, &INFO);

	printf("kalman_evolve (incomplete)\n");

	matrix_free(y);
	matrix_free(A);
	matrix_free(B);

	matrix_free(V_i_c_i);
	matrix_free(V_i_F_i);
	matrix_free(V_i_H_i);
}

void kalman_observe(kalman_t* kalman, matrix_t* G_i, matrix_t* o_i, cov_t* C_i) {

	farray_append(kalman->steps,kalman->current);
	kalman->current = NULL;

	printf("kalman_observe (incomplete)\n");
}


/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

#if false

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
#endif

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

