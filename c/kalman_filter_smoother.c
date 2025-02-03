/*
 * ultimatekalman.c
 *
 * (C) Sivan Toledo, 2022
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef PARALLEL
#include <pthread.h>

//#include "parallel_for_c.h"
void parallel_for_c(void* kalman, void** helper, size_t l, size_t n, size_t block_size, void (*func)(void*, void**, size_t, size_t, size_t));
void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride);
#endif

#ifdef _WIN32
// for "unused" attribute
#define __attribute__(x)
#include <float.h>
// string.h for memcpy
#else
#include <unistd.h>
#endif

/******************************************************************************/
/* BLAS AND LAPACK DECLARATIONS                                               */
/******************************************************************************/

#if defined(BUILD_MEX) && defined(BUILD_MATLAB)
#define blas_int_t mwSignedIndex
#define HAS_BLAS_INT
#else
#define blas_int_t int32_t
#endif

#ifdef BUILD_BLAS_INT
#define blas_int_t BUILD_BLAS_INT
#define HAS_BLAS_INT
#endif

#ifdef BUILD_MKL
#define HAS_BLAS_H
#define HAS_LAPACK_H
#define blas_int_t MKL_INT
#define HAS_BLAS_INT
#ifdef BUILD_BLAS_UNDERSCORE
#undef BUILD_BLAS_UNDERSCORE
#endif
#include <mkl.h>
#endif

#ifndef HAS_BLAS_INT
#define blas_int_t int32_t
#define HAS_BLAS_INT
#endif

#ifdef BUILD_BLAS_H
#define HAS_BLAS_H
#include <blas.h>
#endif

#ifdef BUILD_LAPACK_H
#define HAS_LAPACK_H
#include <lapack.h>
#endif

#ifndef HAS_LAPACK_H
// Microsoft's cl does not support #warning
//#warning "LAPACK subroutines defined in ultimatekalman.c, no header file"
void
#ifdef BUILD_BLAS_UNDERSCORE
     dormqr_
#else
     dormqr
#endif
		(
    char const* side, char const* trans,
		blas_int_t const* m, blas_int_t const* n, blas_int_t const* k,
    double const* A, blas_int_t const* lda,
    double const* tau,
    double* C, blas_int_t const* ldc,
    double* work, blas_int_t const* lwork,
		blas_int_t* info
#ifdef BUILD_BLAS_STRLEN_END
    , size_t, size_t
#endif
);

void
#ifdef BUILD_BLAS_UNDERSCORE
     dgeqrf_
#else
     dgeqrf
#endif
		(
		blas_int_t const* m, blas_int_t const* n,
    double* A, blas_int_t const* lda,
    double* tau,
    double* work, blas_int_t const* lwork,
		blas_int_t* info );

void
#ifdef BUILD_BLAS_UNDERSCORE
     dtrtrs_
#else
     dtrtrs
#endif
    (
    char const* uplo, char const* trans, char const* diag,
		blas_int_t const* n, blas_int_t const* nrhs,
    double const* A, blas_int_t const* lda,
    double* B, blas_int_t const* ldb,
		blas_int_t* info
#ifdef BUILD_BLAS_STRLEN_END
    , size_t, size_t, size_t
#endif
);

#endif

#ifndef HAS_BLAS_H
// Microsoft's cl does not support #warning
//#warning "BLAS subroutines defined in ultimatekalman.c, no header file"
void
#ifdef BUILD_BLAS_UNDERSCORE
     dgemm_
#else
     dgemm
#endif
          (char const* transa, char const* transb,
           blas_int_t const* m, blas_int_t const* n, blas_int_t const* k,
					 double* ALPHA,
					 double* A, blas_int_t const* LDA,
					 double* B, blas_int_t const* LDB,
					 double* BETA,
					 double* C, blas_int_t const* LDC
#ifdef BUILD_BLAS_STRLEN_END
           , size_t, size_t
#endif
					 );
#endif

#define ULTIMATEKALMAN_C
#include "kalman_filter_smoother.h"



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

//static int debug = 0;

/******************************************************************************/
/* UTILITIES                                                                  */
/******************************************************************************/

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

static double NaN = 0.0 / 0.0;

#if defined(BUILD_WIN32_GETTIMEOFDAY) && defined(_WIN32)
#include <windows.h>
typedef SSIZE_T ssize_t;

static
int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (long) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);
    return 0;
}
#else
#include <sys/time.h>
#endif

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
	(A->elements)[ j*(A->ld) + i ] = v;
}

double matrix_get(matrix_t* A, int32_t i, int32_t j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < A->row_dim);
	assert(j < A->col_dim);
	return (A->elements)[ j*(A->ld) + i ];
}
#endif

// this replaces the matrix by its top-left block
void matrix_mutate_chop(matrix_t* A, int32_t rows, int32_t cols) {
	assert(rows >= 0);
	assert(rows >= 0);
	assert(rows <= A->row_dim);
	assert(cols <= A->col_dim);

	A->row_dim = rows;
	A->col_dim = cols;
}

void matrix_mutate_triu(matrix_t* A) {
	int32_t n = matrix_rows(A);
	int32_t i,j;

	for (j=0; j<n; j++) {
		for (i=j+1; i<n; i++) {
			matrix_set(A,i,j,0.0);
		}
	}
}

/*
 * Creates a matrix with undefined elements
 */
matrix_t* matrix_create(int32_t rows, int32_t cols) {
	matrix_t* A = malloc(sizeof(matrix_t));
	assert( A!= NULL );
	A->row_dim = rows;
	A->col_dim = cols;
	A->ld      = rows;
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
int32_t matrix_ld  (matrix_t* A) { return A->ld;      }

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
matrix_t* matrix_create_identity(int32_t rows, int32_t cols) {
	int32_t i,j;

	matrix_t* identity = matrix_create(rows,cols);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(identity,i,j, 0.0);
		}
		matrix_set(identity,i,i,1.0);
	}

	/*
	for (i=0; i<MIN(rows,cols); i++) {
		matrix_set(identity,i,i,1.0);
	}
	*/

	return identity;
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
matrix_t* matrix_create_subtract(matrix_t* A, matrix_t* B) {
	int32_t i,j;

	assert(A != NULL);
	assert(B != NULL);
	assert(matrix_rows(A) == matrix_rows(B));
	assert(matrix_cols(A) == matrix_cols(B));

	matrix_t* C = matrix_create_copy(A);

	for (i=0; i<matrix_rows(A); i++) {
		for (j=0; j<matrix_cols(A); j++) {
			matrix_set(C,i,j, matrix_get(A,i,j) - matrix_get(B,i,j));
		}
	}

	return C;
}
matrix_t* matrix_create_add(matrix_t* A, matrix_t* B) {
	int32_t i,j;

	assert(A != NULL);
	assert(B != NULL);
	assert(matrix_rows(A) == matrix_rows(B));
	assert(matrix_cols(A) == matrix_cols(B));

	matrix_t* C = matrix_create_copy(A);

	for (i=0; i<matrix_rows(A); i++) {
		for (j=0; j<matrix_cols(A); j++) {
			matrix_set(C,i,j, matrix_get(A,i,j) + matrix_get(B,i,j));
		}
	}

	return C;
}

matrix_t* matrix_create_transpose(matrix_t* A) {
	int32_t i,j;

	matrix_t* AT = matrix_create(matrix_cols(A),matrix_rows(A));

	for (i=0; i<matrix_rows(AT); i++) {
		for (j=0; j<matrix_cols(AT); j++) {
			matrix_set(AT,i,j, matrix_get(A,j,i));
		}
	}
	return AT;
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

void matrix_mutate_copy(matrix_t* C, matrix_t* A) {
	if (A==NULL && C==NULL) return;
	assert(A != NULL);
	assert(C != NULL);
	assert(matrix_rows(A) == matrix_rows(C));
	assert(matrix_cols(A) == matrix_cols(C));

	int i,j;

	int32_t rows  = (A->row_dim);
	int32_t cols  = (A->col_dim);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(A,i,j));
		}
	}
}

matrix_t* matrix_create_copy(matrix_t* A) {
	if (A==NULL) return NULL;

	int32_t rows  = (A->row_dim);
	int32_t cols  = (A->col_dim);

	matrix_t* C = matrix_create(rows,cols);

	matrix_mutate_copy(C,A);
	/*
	int32_t i,j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(A,i,j));
		}
	}
  */
	return C;
}


matrix_t* matrix_create_sub(matrix_t* A, int32_t first_row, int32_t rows, int32_t first_col, int32_t cols) {
	int i,j;

	assert(first_row >= 0);
	assert(first_col >= 0);
	assert(first_row+rows <= A->row_dim);
	assert(first_col+cols <= A->col_dim);

	matrix_t* C = matrix_create(rows,cols);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(A,first_row+i,first_col+j));
		}
	}

	return C;
}

void matrix_mutate_copy_sub(matrix_t* C, int32_t first_row, int32_t first_col, matrix_t* A) {
	if (A==NULL && C==NULL) return;
	assert(A != NULL);
	assert(C != NULL);
	assert(matrix_rows(C) >= first_row + matrix_rows(A));
	assert(matrix_cols(C) >= first_col + matrix_cols(A));

	int i,j;

	int32_t rows  = (A->row_dim);
	int32_t cols  = (A->col_dim);

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,first_row+i,first_col+j,matrix_get(A,i,j));
		}
	}
}

matrix_t* matrix_create_vconcat(matrix_t* A, matrix_t* B) {
	int i,j;

	int32_t Arows = (A == NULL ? 0 : A->row_dim);
	int32_t Acols = (A == NULL ? 0 : A->col_dim);
	int32_t Brows = (B == NULL ? 0 : B->row_dim);
	int32_t Bcols = (B == NULL ? 0 : B->col_dim);

	//if (debug) printf("vconcat %d %d %d %d\n",Arows,Acols,Brows,Bcols);

	if (Arows + Brows == 0) return NULL;

	int32_t rows = Arows+Brows;
	int32_t cols = MAX(Acols,Bcols); // one could be zero

	//if (debug) printf("vconcat %d %d %d\n",rows,cols,Arows);
	//assert( A->col_dim == B->col_dim );

	matrix_t* C = matrix_create(rows,cols);

	for (i=0; i<Arows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(A,i,j));
		}
	}

	for (   ; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(C,i,j,matrix_get(B,i-Arows,j));
		}
	}

	//if (debug) printf("vconcat done\n");

	return C;
}

void matrix_mutate_scale(matrix_t* A, double s) {
	int32_t i,j;

	int32_t rows = matrix_rows(A);
	int32_t cols = matrix_cols(A);

	for (j=0; j<cols; j++) {
		for (i=0; i<rows; i++) {
			matrix_set(A,i,j, s * matrix_get(A,i,j));
		}
	}
}
matrix_t* matrix_create_multiply(matrix_t* A, matrix_t* B) {
	assert(A != NULL);
	assert(B != NULL);
	assert(matrix_cols(A) == matrix_rows(B));

	matrix_t* output = matrix_create_constant(matrix_rows(A),matrix_cols(B),0.0);
	matrix_mutate_gemm(1, A, B, 0, output);

	return output;
}

// returns TAU
// nor for flat matrices
matrix_t* matrix_create_mutate_qr(matrix_t* A) {
	assert(A != NULL);
	//int32_t i,j;

	int32_t rows = matrix_rows(A);
	int32_t cols = matrix_cols(A);

	assert(rows >= cols);

	blas_int_t M,N,LDA,LWORK,INFO;
	double     WORK_SCALAR;

	M   = matrix_rows(A);
	N   = matrix_cols(A);
	LDA = matrix_ld  (A);

	matrix_t* TAU = matrix_create(N,1);

	LWORK = -1; // tell lapack to compute the size of the work area required
	//if (debug) printf("dgeqrf: M=%d N=%d LDA=%d LWORK=%d\n",M,N,LDA,LWORK);

#ifdef BUILD_BLAS_UNDERSCORE
  dgeqrf_
#else
  dgeqrf
#endif
	      (&M, &N, A->elements, &LDA, TAU->elements, &WORK_SCALAR, &LWORK, &INFO);
	if (INFO != 0) printf("dgeqrf INFO=%d\n",INFO);
	assert(INFO==0);

	LWORK = (blas_int_t) WORK_SCALAR;
	matrix_t* WORK = matrix_create(LWORK,1);
	//printf("dgeqrf requires %f (%d) words in WORK\n",WORK_SCALAR,LWORK);
#ifdef BUILD_BLAS_UNDERSCORE
  dgeqrf_
#else
  dgeqrf
#endif
	      (&M, &N, A->elements, &LDA, TAU->elements, WORK->elements, &LWORK, &INFO);
	if (INFO != 0) printf("dgeqrf INFO=%d\n",INFO);
	assert(INFO==0);
	matrix_free(WORK);

	return TAU;
}

// mutates C
void matrix_mutate_apply_qt(matrix_t* QR, matrix_t* TAU, matrix_t* C) {
	assert(QR != NULL);

	blas_int_t M,N,K,LDA,LDC,LWORK,INFO;
	double     WORK_SCALAR;

	//int32_t rows = matrix_rows(A);
	//int32_t cols = matrix_cols(A);

	LWORK = -1; // tell lapack to compute the size of the work area required

	M = matrix_rows(C);
	N = matrix_cols(C);
	K = matrix_cols(QR); // number of reflectors
	LDA = matrix_ld(QR);
	LDC = matrix_ld(C);
	//printf("dormqr M=%d N=%d K=%d LDA=%d LDC=%d\n",M,N,K,LDA,LDC);

#ifdef BUILD_BLAS_UNDERSCORE
  dormqr_
#else
  dormqr
#endif
        ("L", "T", &M, &N, &K, QR->elements, &LDA, TAU->elements, C->elements, &LDC, &WORK_SCALAR, &LWORK, &INFO
#ifdef BUILD_BLAS_STRLEN_END
         ,1,1
#endif
			  );
	if (INFO != 0) printf("dormqr INFO=%d\n",INFO);
	assert(INFO==0);

	LWORK = (blas_int_t) WORK_SCALAR;
	matrix_t* WORK = matrix_create(LWORK,1);

	//if (debug) printf("dormqr requires %f (%d) words in WORK\n",WORK_SCALAR,LWORK);

	// left (pre) multiplication by Q^T
#ifdef BUILD_BLAS_UNDERSCORE
  dormqr_
#else
  dormqr
#endif
	      ("L", "T", &M, &N, &K, QR->elements, &LDA, TAU->elements, C->elements, &LDC, WORK->elements, &LWORK, &INFO
#ifdef BUILD_BLAS_STRLEN_END
         ,1,1
#endif
			  );
	if (INFO != 0) printf("dormqr INFO=%d\n",INFO);
	assert(INFO==0);

	matrix_free(WORK);
}

// mutates b
void matrix_mutate_trisolve(matrix_t* U, matrix_t* b) {
	assert(U != NULL);
	assert(b != NULL);
	assert(matrix_rows(U) == matrix_cols(U));
	assert(matrix_rows(b) == matrix_cols(U));

	blas_int_t N,NRHS,LDA,LDB,INFO;

	//if (debug) printf("**** %d %d %d %d\n",matrix_rows(A),matrix_ld(A),matrix_rows(kalman->current->state),matrix_ld(kalman->current->state));

	N = matrix_cols(U);
	LDA = matrix_ld(U);

	NRHS = matrix_cols(b);
	LDB = matrix_ld(b);

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("dtrtrs N=%d NRHS=%d LDA=%d LDB=%d\n",N,NRHS,LDA,LDB);
#endif

	//if (debug) matrix_print(kalman->current->Rdiag,NULL);
	//if (debug) matrix_print(state,NULL);

	// upper triangular, no transpose, diagonal is general (not unit)
#ifdef BUILD_BLAS_UNDERSCORE
     dtrtrs_
#else
     dtrtrs
#endif
           ("U","N","N", &N, &NRHS, U->elements, &LDA, b->elements, &LDB, &INFO
#ifdef BUILD_BLAS_STRLEN_END
    ,1,1,1
#endif
);
	if (INFO != 0) {
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("dormqr INFO=%d\n",INFO);
#endif
	}
	assert(INFO==0);
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("dtrtr done\n");
#endif
}


matrix_t* matrix_create_trisolve(matrix_t* U, matrix_t* b) {
	matrix_t* x = matrix_create_copy(b);
	matrix_mutate_trisolve(U,x);
	return x;
}

matrix_t* matrix_create_mldivide(matrix_t* A, matrix_t* B) {
	assert(A != NULL);
	assert(B != NULL);
	assert(matrix_rows(A) == matrix_cols(A));

	matrix_t* R = matrix_create_copy(A);
	matrix_t* Q = matrix_create_mutate_qr(R);

	B = matrix_create_copy(B);
	matrix_mutate_apply_qt(R,Q,B);

	matrix_mutate_triu(R);
	matrix_t* A_inverse =  matrix_create_trisolve(R,B);

	matrix_free(R);
	matrix_free(Q);
	matrix_free(B);

	return A_inverse;

}	

matrix_t* matrix_create_inverse(matrix_t* A) {
	matrix_t* I = matrix_create_identity(matrix_rows(A),matrix_cols(A));
	matrix_t* A_inverse = matrix_create_mldivide(A,I);
	matrix_free(I);

	return A_inverse;
}

// mutates C
void matrix_mutate_gemm(double ALPHA, matrix_t* A, matrix_t* B, double BETA, matrix_t* C) {
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);
	assert(matrix_rows(B) == matrix_cols(A));
	assert(matrix_rows(C) == matrix_rows(A));
	assert(matrix_cols(C) == matrix_cols(B));

	blas_int_t N,M,K,LDA,LDB,LDC;

	//if (debug) printf("**** %d %d %d %d\n",matrix_rows(A),matrix_ld(A),matrix_rows(kalman->current->state),matrix_ld(kalman->current->state));

	M = matrix_rows(C);
	N = matrix_cols(C);
	K = matrix_cols(A);
	LDA = matrix_ld(A);
	LDB = matrix_ld(B);
	LDC = matrix_ld(C);

	//if (debug) printf("dtrtrs N=%d NRHS=%d LDA=%d LDB=%d\n",N,NRHS,LDA,LDB);

	// no transpose
#ifdef BUILD_BLAS_UNDERSCORE
  dgemm_
#else
  dgemm
#endif
	     ("N","N", &M, &N, &K, &ALPHA, A->elements, &LDA, B->elements, &LDB, &BETA, C->elements, &LDC
#ifdef BUILD_BLAS_STRLEN_END
        ,1,1
#endif
			 );
	//if (debug) printf("dtrtr done\n");
}

#ifdef BUILD_MEX
matrix_t* matrix_create_from_mxarray(const mxArray* mx) {
	if (mx==NULL) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:assertion","matrix_create_from_mxarray: argument must not be null");
	}
	int32_t rows = (int32_t) mxGetM(mx);
	int32_t cols = (int32_t) mxGetN(mx);

	if (rows==0 || cols==0) return NULL;

	//if (debug) printf("== create from mx %d %d\n",rows,cols);

	matrix_t* A = matrix_create(rows,cols);

	memcpy( A->elements, mxGetPr(mx), rows*cols*sizeof(double));

	return A;
}

mxArray* matrix_copy_to_mxarray(matrix_t* A) {
	int32_t i,j;

	assert(A != NULL);

	//printf("to mx array !!!\n");
	//printf("  matrix_t* = %08x\n",A);
	//printf("  rows %d\n",A->row_dim);
	//printf("  cols %d\n",A->col_dim);

	int32_t rows = matrix_rows(A);
	int32_t cols = matrix_cols(A);

	//printf("to mx array %d by %d\n",rows,cols);

	mxArray* mx = mxCreateDoubleMatrix(rows,cols,mxREAL);
	assert(mx != NULL);

	double* out = mxGetPr(mx);
	assert(out != NULL);

	for (j=0; j<cols; j++) {
		for (i=0; i<rows; i++) {
			out[ (j*rows) + i ] = matrix_get(A,i,j);
		}
	}

	//printf("  mxArray* = %08x\n",mx);
	return mx;
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

		// if (debug) printf("append l=%lld p=%lld\n",logical_size,array_size);

		if (logical_size <= array_size/2) { // can just shift back less than half the array
			for (i=0; i<logical_size; i++) {
				(a->elements)[ i ] = (a->elements)[ (a->start) + i ];
			}
			a->start = 0;
			a->end   = logical_size-1;
			//if (debug) printf("farray shifted back, physical size %lld pointers, logical size %lld\n",a->array_size,logical_size);
		} else { // array is currently more than half full, realloc at a larger size
			a->array_size *= 2;
			a->elements = realloc( a->elements, (a->array_size)*sizeof(void*) );
			assert( a->elements != NULL);
			//if (debug) printf("farray reallocated, new physical size is %lld pointers, logical size %lld\n",a->array_size,logical_size);
		}
	}

	assert( a->end < (a->array_size)-1 );

	(a->end)++;
	(a->elements)[ a->end ] = v;
}

static void* farray_drop_first(farray_t* a) {
	void* r = (a->elements)[ a->start ];
	(a->first)++;
	(a->start)++;
	return r;
}

static void* farray_drop_last(farray_t* a) {
	void* r = (a->elements)[ a->end ];
	(a->end)--;
	return r;
}

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

matrix_t* cov_weigh(matrix_t* cov, char cov_type, matrix_t* A) {
	blas_int_t M,N,K,LDA,LDB,LDC,NRHS,INFO;
	double ALPHA, BETA;

	assert(A != NULL);
	assert(cov != NULL);

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("cov(%c) ",cov_type);
	matrix_print(cov,"%.3e");
	printf("A ");
	matrix_print(A,"%.3e");
#endif

	matrix_t* WA = NULL;

	int32_t i,j, rows, cols;

	switch (cov_type) {
	case 'W':
		//if (debug) printf("cov W %d %d %d %d\n",matrix_cols(cov),matrix_rows(cov),matrix_cols(A),matrix_rows(A));

		assert(matrix_cols(cov) == matrix_rows(A));

		WA = matrix_create_constant(matrix_rows(A), matrix_cols(A),0.0);

		M = matrix_rows(WA);
		N = matrix_cols(WA);
		K = matrix_cols(cov);
		LDA = matrix_ld(cov);
		LDB = matrix_ld(A);
		LDC = matrix_ld(WA);
		ALPHA = 1.0;
		BETA  = 0.0;
#ifdef BUILD_BLAS_UNDERSCORE
    dgemm_
#else
    dgemm
#endif
		     ("n","N",&M,&N,&K,&ALPHA, cov->elements, &LDA, A->elements, &LDB, &BETA, WA->elements, &LDC
#ifdef BUILD_BLAS_STRLEN_END
          ,1,1
#endif
			 );
		break;
	case 'U': // cov and an upper triangular matrix that we need to solve with
	case 'F': // same representation; in Matlab, we started from explicit cov and factored

		//if (debug) printf("cov U %d %d %d %d\n",matrix_cols(cov),matrix_rows(cov),matrix_cols(A),matrix_rows(A));

		// is this correct for 'F'?
		WA = matrix_create_trisolve(cov,A);
		/*
		assert(matrix_cols(cov) == matrix_rows(A));

		WA = matrix_create_copy(A);

		N = matrix_cols(cov);
		LDA = matrix_ld(cov);

		NRHS = matrix_cols(WA);
		LDB = matrix_ld(WA);

		//if (debug) printf("dtrtrs N=%d NRHS=%d LDA=%d LDB=%d\n",N,NRHS,LDA,LDB);
		//if (debug) printf("dtrtrs A %08x state %08x\n",step->Rdiag->elements, state->elements);

		//if (debug) matrix_print(step->Rdiag,NULL);
		//if (debug) matrix_print(kalman->current->state,NULL);

		// upper triangular, no transpose, diagonal is general (not unit)
		dtrtrs("U","N","N", &N, &NRHS, cov->elements, &LDA, WA->elements, &LDB, &INFO);
		*/
		break;
	case 'w':
		assert(matrix_rows(cov) == matrix_rows(A));

		rows = matrix_rows(A);
		cols = matrix_cols(A);

		WA = matrix_create_constant(rows, cols,0.0);

		for (j=0; j<cols; j++) {
			for (i=0; i<rows; i++) {
				matrix_set(WA,i,j, matrix_get(cov,i,0) * matrix_get(A,i,j) );
			}
		}
		break;
	default:
		printf("unknown covariance-matrix type %c\n",cov_type);
		assert( 0 );
		WA = matrix_create_constant(matrix_rows(A), matrix_cols(A), NaN);
		break;
	}

	//if (debug) printf("WA ");
	//if (debug) matrix_print(WA,"%.3e");

	return WA;
}
// SUPPORT 'W','C' only
matrix_t* explicit(matrix_t* cov, char type) {
	assert(type == 'W' || type == 'C' || type == 'U' || type == 'F');
	if (type == 'W') {
		matrix_t* WT = matrix_create_transpose(cov);
		matrix_t* WTW  = matrix_create_multiply(WT,cov);
		matrix_t* C = matrix_create_inverse(WTW);
		matrix_free(WT);
		matrix_free(WTW);
		return C;
	}
	if (type == 'U' || type == 'F') {
		//printf("L = ");
		//matrix_print(cov,"%.3e");
		matrix_t* LT = matrix_create_transpose(cov);
		matrix_t* C  = matrix_create_multiply(cov,LT);
		matrix_free(LT);
		//printf("C = ");
		//matrix_print(C,"%.3e");
		return C;
	}
	if (type == 'C') {
		matrix_t* C = matrix_create_copy(cov);
		return C;
	}
}

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

static step_t* step_create() {
	step_t* s = malloc(sizeof(step_t));
	s->step      = -1;
	s->dimension = -1;
	
	s->C_type    = 0;
	s->K_type    = 0;


	//s->H = NULL;
	s->F = NULL;
	
	s->predictedState = NULL;
	s->predictedCovariance = NULL;
	s->assimilatedState = NULL;
	s->assimilatedCovariance = NULL;
	s->smoothedState = NULL;
	s->smoothedCovariance = NULL;
	s->state = NULL;
	s->covariance = NULL;

	assert( s != NULL );
	return s;
}

void step_free(step_t* s) {
	
	//matrix_free(s->H);
	matrix_free(s->F);
	
	matrix_free(s->predictedState);
	matrix_free(s->predictedCovariance);
	matrix_free(s->assimilatedState);
	matrix_free(s->assimilatedCovariance);
	matrix_free(s->smoothedState);
	matrix_free(s->smoothedCovariance);

	// state and covariance are aliases in this implementation

	free(s);
}

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

kalman_t* kalman_create() {
#ifdef PARALLEL
	char* nthreads_string = getenv("NTHREADS");
	int nthreads = 0;
	if (nthreads_string != NULL && sscanf(nthreads_string,"%d",&nthreads)==1) {
		kalman_parallel_init(nthreads);
		printf("limiting to %d threads/cores\n",nthreads);
	}
#endif

	kalman_t* kalman = malloc(sizeof(kalman_t));
	assert( kalman != NULL );
	kalman->steps   = farray_create();
	kalman->current = NULL;
	return kalman;
}

void kalman_free(kalman_t* kalman) {
	//printf("waning: kalman_free not fully implemented yet (steps not processed)\n");

	while (farray_size(kalman->steps) > 0) {
		step_t* i = farray_drop_last(kalman->steps);
		step_free(i);
	}

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

void kalman_evolve(kalman_t* kalman, int32_t n_i, matrix_t* H_i, matrix_t* F_i, matrix_t* c_i, matrix_t* K_i, char K_type) {
	kalman->current = step_create();
	kalman->current->dimension = n_i;

	if (farray_size(kalman->steps)==0) {
		//if (debug) printf("kalman_evolve first step\n");
		kalman->current->step = 0;
		return;
	}

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_evolve n_i = %d\n",n_i);
	printf("kalman_evolve H_i = %08x %d %d\n",(unsigned int) H_i,(int) matrix_rows(H_i),(int) matrix_cols(H_i));
	printf("kalman_evolve F_i = %08x %d %d\n",(unsigned int) F_i,(int) matrix_rows(F_i),(int) matrix_cols(F_i));
	printf("kalman_evolve i_i = %08x %d %d\n",(unsigned int) c_i,(int) matrix_rows(c_i),(int) matrix_cols(c_i));
	printf("kalman_evolve K_i = %08x %d %d\n",(unsigned int) K_i,(int) matrix_rows(K_i),(int) matrix_cols(K_i));
	printf("cov type %c\n",K_type);
#endif

	step_t* imo = farray_get_last( kalman->steps );
	assert(imo != NULL);
	kalman->current->step = (imo->step) + 1;

	//if (debug) printf("kalman_evolve step = %d\n",kalman->current->step);

	assert(H_i!=NULL);
	assert(F_i!=NULL);
	assert(c_i!=NULL);

	// we assume H_i is an identity, need to check in the final code
	//kalman->current->H	= matrix_create_copy(H_i);
	kalman->current->F 	= matrix_create_copy(F_i);

	//printf("imo->assimilatedState = ");
	//matrix_print(imo->assimilatedState,"%.3e");
	matrix_t* predictedState = matrix_create_multiply( F_i, imo->assimilatedState );

	// with H_i
	//matrix_t* Hinv = matrix_create_inverse(H_i);
	//matrix_t* HinvTrans = matrix_create_transpose( Hinv );
	//matrix_t* predictedPlusCi = matrix_create_add(predictedState,c_i);
	//matrix_free(predictedState);
	//kalman->current->predictedState = matrix_create_multiply(Hinv, predictedPlusCi);
	//matrix_free(predictedPlusCi);

	//matrix_t* K_i_explicit = explicit(K_i,K_type);
	//matrix_t* t1 = matrix_create_multiply(Hinv,K_i_explicit);
	//matrix_t* t2 = matrix_create_multiply(t1,HinvTrans);

	//matrix_t* t3 = matrix_create_multiply( Hinv, F_i );
	//matrix_t* t4 = matrix_create_multiply( t3, imo->assimilatedCovariance );
	//matrix_t* F_iTrans = matrix_create_transpose( F_i );
	//matrix_t* t5 = matrix_create_multiply( t4, F_iTrans );
	//matrix_t* t6 = matrix_create_multiply( t5, HinvTrans );

	//matrix_free( t6 );
	//matrix_free( F_iTrans );
	//matrix_free( t5 );
	//matrix_free( t4 );
	//matrix_free( t3 );
	//matrix_free( t2 );
	//matrix_free( t1 );
	//matrix_free( HinvTrans );
	//matrix_free( Hinv );

	kalman->current->predictedState = matrix_create_add(predictedState,c_i);
	matrix_free(predictedState);

	matrix_t* K_i_explicit = explicit(K_i,K_type);

	matrix_t* t4 = matrix_create_multiply( F_i, imo->assimilatedCovariance );
	matrix_t* F_iTrans = matrix_create_transpose( F_i );
	matrix_t* t5 = matrix_create_multiply( t4, F_iTrans );

	kalman->current->predictedCovariance = matrix_create_add( t5, K_i_explicit );

	matrix_free( F_iTrans );
	matrix_free( t5 );
	matrix_free( t4 );

	kalman->current->state      = kalman->current->predictedState      ;
	kalman->current->covariance = kalman->current->predictedCovariance ;
}

void kalman_observe(kalman_t* kalman, matrix_t* G_i, matrix_t* o_i, matrix_t* C_i, char C_type) {

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe!\n");
#endif

	assert(kalman != NULL);
	assert(kalman->current != NULL);
	int32_t n_i = kalman->current->dimension;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe %d\n",(int) kalman->current->step);
#endif

	// matrix_t* W_i_G_i = NULL;
	// matrix_t* W_i_o_i = NULL;

	if (kalman->current->step == 0) {
		//printf("filter_smoother step 0 observation cov-type %c\n",C_type);
		matrix_t* W_i_G_i = cov_weigh(C_i,C_type,G_i);
		matrix_t* W_i_o_i = cov_weigh(C_i,C_type,o_i);

		//matrix_print(W_i_G_i,"%.3e");
		//matrix_print(W_i_o_i,"%.3e");

		matrix_t* tau = matrix_create_mutate_qr( W_i_G_i );
		matrix_mutate_apply_qt( W_i_G_i, tau, W_i_o_i );

		matrix_mutate_chop(W_i_G_i, n_i, n_i); // it might have been tall
		matrix_mutate_chop(W_i_o_i, n_i, 1  ); // it might have been tall
		matrix_mutate_triu(W_i_G_i);

		kalman->current->assimilatedState =  matrix_create_trisolve(W_i_G_i,W_i_o_i);

		matrix_t* R_trans = matrix_create_transpose(W_i_G_i);
		matrix_t* R_trans_R = matrix_create_multiply( R_trans, W_i_G_i );

		kalman->current->assimilatedCovariance = matrix_create_inverse(R_trans_R);

		//matrix_print(kalman->current->assimilatedCovariance,"%.3e");
		//matrix_print(kalman->current->assimilatedState,"%.3e");

		matrix_free( R_trans_R );
		matrix_free( R_trans );
		matrix_free( tau );
		matrix_free( W_i_G_i );
		matrix_free( W_i_o_i );

		kalman->current->state      = kalman->current->assimilatedState      ;
		kalman->current->covariance = kalman->current->assimilatedCovariance ;

		farray_append(kalman->steps,kalman->current);
		return;
	}

	if (o_i == NULL) {
		kalman->current->assimilatedState      = matrix_create_copy( kalman->current->predictedState      );
		kalman->current->assimilatedCovariance = matrix_create_copy( kalman->current->predictedCovariance );

	} else {
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("C_i(%c) ",C_type);
		matrix_print(C_i,"%.3e");
		printf("G_i ");
		matrix_print(G_i,"%.3e");
		printf("o_i ");
		matrix_print(o_i,"%.3e");
#endif
		// W_i_G_i = cov_weigh(C_i,C_type,G_i);
		// W_i_o_i = cov_weigh(C_i,C_type,o_i);

		matrix_t* predictedObservations = matrix_create_multiply( G_i, kalman->current->predictedState );

		matrix_t* G_i_trans = matrix_create_transpose( G_i );
		matrix_t* t1 = matrix_create_multiply( G_i, kalman->current->predictedCovariance );
		matrix_t* t2 = matrix_create_multiply( t1, G_i_trans );
		matrix_t* C_i_explicit = explicit( C_i, C_type );
		matrix_t* S = matrix_create_add( t2, C_i_explicit );

		matrix_t* t4 = matrix_create_multiply( kalman->current->predictedCovariance, G_i_trans );
		matrix_t* S_inv = matrix_create_inverse( S );
		matrix_t* gain = matrix_create_multiply( t4, S_inv );

		matrix_t* innovation = matrix_create_subtract( o_i, predictedObservations );

		matrix_t* gain_innovation = matrix_create_multiply( gain, innovation );

		kalman->current->assimilatedState = matrix_create_add( kalman->current->predictedState, gain_innovation );

		matrix_t* t5 = matrix_create_multiply( gain, G_i );
		matrix_t* t6 = matrix_create_multiply( t5, kalman->current->predictedCovariance );

		kalman->current->assimilatedCovariance = matrix_create_subtract( kalman->current->predictedCovariance, t6 );

		matrix_free( t6 );
		matrix_free( gain_innovation );
		matrix_free( innovation );
		matrix_free( gain );
		matrix_free( S_inv );
		matrix_free( t5 );
		matrix_free( t4 );
		matrix_free( S );
		matrix_free( C_i_explicit );
		//matrix_free( t3 );
		matrix_free( t2 );
		matrix_free( t1 );
		matrix_free( G_i_trans );
		matrix_free( predictedObservations );

	}

	kalman->current->state      = kalman->current->assimilatedState      ;
	kalman->current->covariance = kalman->current->assimilatedCovariance ;

	farray_append(kalman->steps,kalman->current);

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_observe done\n");
#endif
}

matrix_t* kalman_estimate(kalman_t* kalman, int64_t si) {
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_estimate\n");
#endif

	if (farray_size(kalman->steps) == 0) return NULL;

	if (si < 0) si = farray_last_index(kalman->steps);
	step_t* step = farray_get(kalman->steps,si);

	if (step->state == NULL) return matrix_create_constant(step->dimension,1,NaN);

	return matrix_create_copy(step->state);
}

void kalman_smooth(kalman_t* kalman) {
	smooth(kalman);
}

void smooth(kalman_t* kalman) {
	if (farray_size(kalman->steps) == 0) return;

	int64_t si;
	int64_t last  = farray_last_index (kalman->steps);
	int64_t first = farray_first_index(kalman->steps);

	step_t* i = farray_get(kalman->steps,last);

	i->smoothedState      = matrix_create_copy( i->assimilatedState      );
	i->smoothedCovariance = matrix_create_copy( i->assimilatedCovariance );

	i->state      = i->smoothedState      ;
	i->covariance = i->smoothedCovariance ;

	printf("smooth first:last = %d:%d\n",first,last);
	step_t* ipo = i;
	for (si=last-1; si>=first; si--) {
		i = farray_get(kalman->steps,si);

		matrix_t* nextPredictedEstimate    = ipo->predictedState;
		matrix_t* nextPredictedCovariance  = ipo->predictedCovariance;

		matrix_t* nextPredictedCovarianceInv  = matrix_create_inverse( nextPredictedCovariance );

		matrix_t* nextSmoothEstimate       = ipo->smoothedState;
		matrix_t* nextSmoothCovariance     = ipo->smoothedCovariance;

		//matrix_t* nextEvolutionMatrix      = matrix_create_mldivide( ipo->H, ipo->F );
		matrix_t* nextEvolutionMatrix      = ipo->F;
		matrix_t* nextEvolutionMatrixTrans = matrix_create_transpose( nextEvolutionMatrix );

		matrix_t* assimilatedState         = i->assimilatedState;
		matrix_t* assimilatedCovariance    = i->assimilatedCovariance;

		matrix_t* t1                       = matrix_create_multiply( assimilatedCovariance, nextEvolutionMatrixTrans );
		matrix_t* backwardInnovation       = matrix_create_multiply( t1, nextPredictedCovarianceInv );
		matrix_t* backwardInnovationTrans  = matrix_create_transpose( backwardInnovation );
		matrix_t* t2                       = matrix_create_subtract( nextSmoothEstimate, nextPredictedEstimate );
		matrix_t* t3                       = matrix_create_multiply( backwardInnovation, t2 );
		i->smoothedState                   = matrix_create_add( assimilatedState, t3 );

		matrix_t* t4                       = matrix_create_subtract( nextSmoothCovariance, nextPredictedCovariance );
		matrix_t* t5                       = matrix_create_multiply( backwardInnovation, t4 );
		matrix_t* t6                       = matrix_create_multiply( t5, backwardInnovationTrans );
		i->smoothedCovariance              = matrix_create_add( assimilatedCovariance, t6 );

		matrix_free( t6 );
		matrix_free( t5 );
		matrix_free( t4 );
		matrix_free( t3 );
		matrix_free( t2 );
		matrix_free( t1 );
		//matrix_free( nextEvolutionMatrix );
		matrix_free( nextEvolutionMatrixTrans );
		matrix_free( backwardInnovationTrans );
		matrix_free( backwardInnovation );
		matrix_free( nextPredictedCovarianceInv );

		i->state      = i->smoothedState      ;
		i->covariance = i->smoothedCovariance ;

		ipo = i;
	}

}

char kalman_covariance_type(kalman_t* kalman, int64_t si) { return 'C'; }

matrix_t* kalman_covariance(kalman_t* kalman, int64_t si) {
	
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_covariance\n");
#endif

	if (farray_size(kalman->steps) == 0) return NULL;

	if (si < 0) si = farray_last_index(kalman->steps);
	step_t* step = farray_get(kalman->steps,si);

	int32_t n_i = step->dimension;

	matrix_t* cov = NULL;

	if (step->covariance) {
		cov = matrix_create_copy(step->covariance);
	} else {
		cov = matrix_create_constant(n_i,n_i,NaN);
	}

	return cov;
}


void kalman_forget(kalman_t* kalman, int64_t si) {
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d\n",(int) si);
#endif

	if (farray_size(kalman->steps) == 0) return;

	if (si < 0) si = farray_last_index(kalman->steps) - 1;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d now %d to %d\n",(int) si,(int) farray_first_index(kalman->steps),(int) farray_last_index(kalman->steps));
#endif

	if (si > farray_last_index (kalman->steps) - 1) return; // don't delete last step
	if (si < farray_first_index(kalman->steps)    ) return; // nothing to delete


	while (farray_first_index(kalman->steps) <= si) {
		step_t* step = farray_drop_first(kalman->steps);
		step_free(step);
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("forget new first %d\n",(int) farray_first_index(kalman->steps));
#endif
	}
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d done\n",(int) si);
#endif
}

void kalman_rollback(kalman_t* kalman, int64_t si) {
	if (farray_size(kalman->steps) == 0) return;

	if (si > farray_last_index (kalman->steps)) return; // we can roll  back even the last step (its observation)
	if (si < farray_first_index(kalman->steps)) return;

	step_t* step;
	do {
		step = farray_drop_last(kalman->steps);
		if (step->step == si) {

			matrix_free(step->smoothedState);
			matrix_free(step->smoothedCovariance);
			matrix_free(step->assimilatedState);
			matrix_free(step->assimilatedCovariance);
			// fix aliases
			step->state      = step->predictedState;
			step->covariance = step->predictedCovariance;

			kalman->current = step;

		} else {
			//printf("rollback dropped step %d\n",step->step);
		}
	} while (step->step > si);
	//printf("rollback to %d new latest %d\n",si,kalman_latest(kalman));
}

static struct timeval begin, end;

matrix_t* kalman_perftest(kalman_t* kalman,
		                      matrix_t* H, matrix_t* F, matrix_t* c, matrix_t* K, char K_type,
		                      matrix_t* G, matrix_t* o,              matrix_t* C, char C_type,
									        int32_t count, int32_t decimation) {

	matrix_t* t = matrix_create_constant(count/decimation, 1, NaN);
	int32_t i,j,n;

	//printf("perftest count %d decimation %d rows %d\n",count,decimation,matrix_rows(t));

	//struct timeval begin, end;
	gettimeofday(&begin, 0);

	j = 0;
	n = matrix_cols(G);

	for (i=0; i<count; i++) {
		//printf("perftest iter %d (j=%d)\n",i,j);
		//if (debug) printf("perftest iter %d (j=%d)\n",i,j);
		kalman_evolve(kalman,n,H,F,c,K,K_type);
		kalman_observe(kalman,G,o,C,C_type);
		matrix_t* e = kalman_estimate(kalman,-1);
		matrix_free(e);
		kalman_forget(kalman,-1);

		if ((i % decimation) == decimation-1) {
			gettimeofday(&end, 0);
			long seconds      = end.tv_sec  - begin.tv_sec;
			long microseconds = end.tv_usec - begin.tv_usec;
			double elapsed    = seconds + microseconds*1e-6;

			assert(j < matrix_rows(t));

			//printf("perftest store i=%d (j=%d) %.3e\n",i,j,elapsed/decimation);
			matrix_set(t,j,0,elapsed / decimation); // average per step
			j = j+1;

			gettimeofday(&begin, 0);
			//begin = end;
		}
	}

	//printf("perftest done\n");
	//printf("  matrix_t* t = %08x\n",t);

	//printf("perftest count %d decimation %d rows %d cols %d\n",count,decimation,matrix_rows(t),matrix_cols(t));

	//matrix_print(t,"%.3e");

	return t;
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
