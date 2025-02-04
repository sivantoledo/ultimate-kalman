/*
 * ultimatekalman.c
 *
 * (C) Sivan Toledo, 2022
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// for getenv
#include <stdlib.h>

#include <math.h>

//#include "parallel_for_c.h"
void parallel_for_c(void* kalman, void* indices, int length, size_t n, size_t block_size, void (*func)(void*, void*, int, size_t, size_t));
int kalman_parallel_init(int number_of_threads);

#ifdef PARALLEL
#define malloc(x) aligned_alloc(64,(x))
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
//#include "ultimatekalman_oddeven.h"

typedef void kalman_t; 

typedef struct matrix_st {
	int32_t row_dim;
	int32_t col_dim;
	int32_t ld;      // leading dimension
	double* elements;
}
kalman_matrix_t
#ifdef ULTIMATEKALMAN_C
,matrix_t
#endif
;

#ifdef NDEBUG
#define matrix_set(A,i,j,v) ( ((A)->elements)[ (j)*((A)->ld) + (i)  ] = (v) )
#define matrix_get(A,i,j)   ( ((A)->elements)[ (j)*((A)->ld) + (i)  ] )
#else
void   matrix_set(kalman_matrix_t* A, int32_t i, int32_t j, double v);
double matrix_get(kalman_matrix_t* A, int32_t i, int32_t j);
#endif
matrix_t* matrix_create_copy(matrix_t* A);


typedef struct step_st {
	kalman_matrix_t* A;
	kalman_matrix_t* TAU;

	kalman_matrix_t* X3; // pad to cache lines
	kalman_matrix_t* X4;
	kalman_matrix_t* X5;
	kalman_matrix_t* X6;
	kalman_matrix_t* X7;
	kalman_matrix_t* X8;
} step_t;


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

int32_t matrix_rows(matrix_t* A);

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
matrix_t* matrix_mutate_subtract(matrix_t* A, matrix_t* B) {
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

void matrix_print(matrix_t* A, char* format) {
	int32_t i,j;

	// TODO print null if A is null?

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
	// printf("dormqr M=%d N=%d K=%d LDA=%d LDC=%d\n",M,N,K,LDA,LDC);

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

/**************************************************************************************************/

// this function is from ChatGPT.
double generateGaussian(double mean, double stddev) {
    // Generate two uniform random numbers in the range (0, 1)
    double u1 = rand() / (RAND_MAX + 1.0);
    double u2 = rand() / (RAND_MAX + 1.0);

    // Apply Box-Muller transform
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

    // Adjust for desired mean and standard deviation
    return z0 * stddev + mean;
}

//kalman_matrix_t* matrix_create_mutate_qr(kalman_matrix_t* A);
//void matrix_mutate_apply_qt(kalman_matrix_t* QR, kalman_matrix_t* TAU, kalman_matrix_t* C);

kalman_matrix_t* generateRandomOrthonormal(int rows, int cols) {
	int i,j;
	int n;

	if (rows==1 || cols == 1) {
		kalman_matrix_t* QR = matrix_create(rows,cols);
		for (i=0; i<rows; i++) {
			for (j=0; j<cols; j++) {
				matrix_set(QR,i,j,generateGaussian(0.0,1.0));
			}
		}
		return QR;
	}

	// we assume that rows >= cols

	kalman_matrix_t* QR = matrix_create(rows,rows);
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			matrix_set(QR,i,j,generateGaussian(0.0,1.0));
		}
	}

	kalman_matrix_t* A = matrix_create_identity(rows,cols);
	kalman_matrix_t* TAU = matrix_create_mutate_qr(QR);
	matrix_mutate_apply_qt(QR,TAU,A);

	matrix_free(TAU);
	matrix_free(QR);

	return A;
}

double times[16];
matrix_t* A;

//void ep_alloc_struct(kalman_t* kalman, step_t* *indices, int n, size_t start, size_t end) {
void ep_alloc_struct(void* kalman_v, void* indices_v, int n, size_t start, size_t end) {
  kalman_t* kalman  = (kalman_t*) kalman_v;
  step_t**  indices = (step_t**)  indices_v;
  
  for (int i = start; i < end; ++i) {
    indices[i] = (step_t*) malloc(sizeof(step_t));
  }
}

void ep_alloc_matrix(void* kalman_v, void* indices_v, int n, size_t start, size_t end) {
  kalman_t* kalman  = (kalman_t*) kalman_v;
  step_t**  indices = (step_t**)  indices_v;

	for (int i = start; i < end; ++i) {
	  indices[i] -> A = matrix_create(2*n,n);
	}
}

void ep_fill_matrix(void*  kalman_v, void* indices_v, int n, size_t start, size_t end) {
  kalman_t* kalman  = (kalman_t*) kalman_v;
  step_t**  indices = (step_t**)  indices_v;

	for (int i = start; i < end; ++i) {
	  for (int r=0; r<2*n; r++) {
	    for (int c=0; c<n; c++) {
	      matrix_set(indices[i] -> A,r,c,(double) (r+c));
	    }
	  }
	}
}

void ep_create(void*  kalman_v, void* indices_v, int n, size_t start, size_t end) {
  kalman_t* kalman  = (kalman_t*) kalman_v;
  step_t**  indices = (step_t**)  indices_v;

	for (int i = start; i < end; ++i) {
	  indices[i] = (step_t*) malloc(sizeof(step_t));
	  //indices[i] -> A = matrix_create_copy(A);
	  indices[i] -> A = matrix_create(2*n,n);
	  for (int r=0; r<2*n; r++) {
	    for (int c=0; c<n; c++) {
	      matrix_set(indices[i] -> A,r,c,(double) (r+c));
	    }
	  }
	}
}

void ep_factor(void* kalman_v, void* indices_v, int n, size_t start, size_t end) {
  kalman_t* kalman  = (kalman_t*) kalman_v;
  step_t**  indices = (step_t**)  indices_v;

	for (int i = start; i < end; ++i) {
		indices[i]->TAU = matrix_create_mutate_qr(indices[i]->A);
	}
}

int main(int argc, char* argv[]) {

	printf("embarrassingly parallel testing starting (start with two arguments, dimension and count)\n");

#ifdef PARALLEL
	char* nthreads_string = getenv("NTHREADS");
	int nthreads = 0;
	if (nthreads_string != NULL && sscanf(nthreads_string,"%d",&nthreads)==1) {
		kalman_parallel_init(nthreads);
		printf("limiting to %d threads/cores\n",nthreads);
	}

	//char* blocksize_string = getenv("TBB_BLOCKSIZE");
	//int blocksize = 0;
	//if (blocksize_string != NULL && sscanf(blocksize_string,"%d",&blocksize)==1) {
	//	kalman_parallel_blocksize(blocksize);
	//	printf("setting blocksize to %d\n",blocksize);
	//}
#endif

	if (argc < 3) {
		printf("usage: performance dimension count\n");
		printf("       using defaults\n");
	}

	int n = 6;
	int k = 100000;

	if (argc >= 2) sscanf(argv[1],"%d",&n);
	if (argc >= 3) sscanf(argv[2],"%d",&k);

	printf("embarrassingly parallel testing smooth dimension=%d step count=%d, starting\n",n,k);

	double t = 0.0;

	struct timeval begin, end;
	long seconds, microseconds;

	A = matrix_create(2*n,n);
	for (int r=0; r<2*n; r++) {
	  for (int c=0; c<n; c++) {
	    matrix_set(A,r,c,generateGaussian(0.0,1.0));
	  }
	}

	gettimeofday(&begin, 0);
	step_t* *indices = (step_t**)malloc(k * sizeof(step_t*));

	//#ifdef PARALLEL
	//parallel_for_c(NULL, indices, n, k, 8, ep_create);
	//#else
	//ep_create(NULL, indices, n, 0, k);
	//#endif

#ifdef PARALLEL
	parallel_for_c(NULL, indices, n, k, 8, ep_alloc_struct);
#else
	ep_alloc_struct(NULL, indices, n, 0, k);
#endif
	
	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[0]          = seconds + microseconds*1e-6;

#ifdef PARALLEL
	parallel_for_c(NULL, indices, n, k, 8, ep_alloc_matrix);
#else
	ep_alloc_matrix(NULL, indices, n, 0, k);
#endif
	
	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[1]          = seconds + microseconds*1e-6;

#ifdef PARALLEL
	parallel_for_c(NULL, indices, n, k, 8, ep_fill_matrix);
#else
	ep_fill_matrix(NULL, indices, n, 0, k);
#endif
	
	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[2]          = seconds + microseconds*1e-6;
	
#ifdef PARALLEL
	parallel_for_c(NULL, indices, n, k, 8, ep_factor);
#else
	ep_factor(NULL, indices, n, 0, k);
#endif

	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[3]          = seconds + microseconds*1e-6;

	/*
	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[2]          = seconds + microseconds*1e-6;

	gettimeofday(&end, 0);
	seconds      = end.tv_sec  - begin.tv_sec;
	microseconds = end.tv_usec - begin.tv_usec;
	times[3]          = seconds + microseconds*1e-6;

	printf("embarrassingly parallel testing took %.2e seconds\n",times[1]);

	printf("performance testing breakdown %.2e %.2e %.2e %.2e (create, factor, nothing, nothing)\n",
			times[0],
			times[1]-times[0],
			times[2]-times[1],
			times[3]-times[2]);
	*/

	printf("embarrassingly parallel testing took %.2e seconds\n",times[1]);

	printf("performance testing breakdown %.2e %.2e %.2e %.2e (alloc_struct, alloc_matrix, fill_matrix, factor)\n",
			times[0],
			times[1]-times[0],
			times[2]-times[1],
			times[3]-times[2]);

	printf("embarrassingly parallel testing done\n");

	return 0;
}


/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
