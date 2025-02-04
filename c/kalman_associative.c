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
int kalman_parallel_init(int number_of_threads);
void parallel_for_c(void* kalman, void** helper, size_t l, size_t n, size_t block_size, void (*func)(void*, void**, size_t, size_t, size_t));
void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride);
#endif

#ifdef PARALLEL
#ifdef MACOS
void* local_aligned_alloc(size_t alignment, size_t size) {
	void* p;
	int result = posix_memalign(&p, alignment, size);
	if (result != 0) return NULL;
	return p;
} 
#define malloc(x) local_aligned_alloc(64,(x))
#else
#define malloc(x) aligned_alloc(64,(x))
#endif
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
#include "kalman_associative.h"



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
	return NULL;
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


	s->H = NULL;
	s->F = NULL;
	s->K = NULL;
	s->c = NULL;

	s->G = NULL;
	s->o = NULL;
	s->C = NULL;

	s->Z = NULL;

	s->A = NULL;
	s->b = NULL;

	s->e = NULL;
	s->J = NULL;

	s->E = NULL;
	s->g = NULL;
	s->L = NULL;
	
	s->state = NULL;
	s->covariance = NULL;

	assert( s != NULL );
	return s;
}

void step_free(step_t* s) {
	
	matrix_free(s->H);
	matrix_free(s->F);
	matrix_free(s->K);
	matrix_free(s->c);

	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->Z);

	matrix_free(s->A);
	matrix_free(s->b);

	matrix_free(s->e);
	matrix_free(s->J);

	matrix_free(s->E);
	matrix_free(s->g);
	matrix_free(s->L);
	
	matrix_free(s->state);
	matrix_free(s->covariance);

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
	//assert(K_i!=NULL);


	// matrix_t* V_i_H_i = cov_weigh(K_i,K_type,H_i);
	// matrix_t* V_i_F_i = cov_weigh(K_i,K_type,F_i);
	// matrix_t* V_i_c_i = cov_weigh(K_i,K_type,c_i);

	// matrix_mutate_scale(V_i_F_i,-1.0);

	kalman->current->H	= matrix_create_copy(H_i);
	kalman->current->F 	= matrix_create_copy(F_i);
	kalman->current->c  = matrix_create_copy(c_i);
	kalman->current->K  = matrix_create_copy(K_i);
	kalman->current->K_type = K_type;
	
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

	if (o_i != NULL) {
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

	kalman->current->G  = matrix_create_copy(G_i);
	kalman->current->o  = matrix_create_copy(o_i);
	kalman->current->C  = matrix_create_copy(C_i);
	kalman->current->C_type = C_type;

#ifdef BUILD_DEBUG_PRINTOUTS
		printf("W_i_G_i ");
		matrix_print(W_i_G_i,"%.3e");
		printf("W_i_o_i ");
		matrix_print(W_i_o_i,"%.3e");
#endif
	}

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

#ifdef PARALLEL
/******************************************************************************/
/* Lock Arrays                                                                */
/******************************************************************************/

#define COLUMNS 10
//#define BLOCKSIZE 1000
#define BLOCKSIZE 10

void parallelInit(void* la_v, void* *helper, size_t length, size_t start, size_t end){
    LockedArray_t* la = (LockedArray_t*) la_v;
	for (int i = start; i < end; i++) {
        la->arrays[i] = (step_t**)malloc(COLUMNS * sizeof(step_t*));
        for (int j = 0; j < COLUMNS; j++) {
            la->arrays[i][j] = (step_t*)NULL;
        }
#ifdef PARALLEL
        pthread_mutex_init(&la->locks[i], NULL);  // Initialize the lock for the row
#endif
    }
}

LockedArray_t* initLockedArray(int k) {
    LockedArray_t* la = (LockedArray_t*)malloc(sizeof(LockedArray_t));
    la->rows = k;
    la->columns = COLUMNS;
    la->arrays = (step_t***)malloc(k * sizeof(step_t**));
#ifdef PARALLEL
    la->locks = (pthread_mutex_t*)malloc(k * sizeof(pthread_mutex_t));
#endif

#ifdef PARALLEL
	parallel_for_c(la, NULL, 0, k, BLOCKSIZE, parallelInit);
#else
	fprintf(stderr,"1\n");
	parallelInit(la, NULL, 0, 0, k); // k was l
	fprintf(stderr,"2\n");
#endif

    return la;
}

int findEmptyColumn(step_t** row, int columns) {
    int attempts = 0;
    while (attempts < columns) {
        int col = rand() % columns;
        if (row[col] == (step_t*)NULL) {
            return col;
        }
        attempts++;
    }
	assert(0);
    return -1;  // No empty column found, need to increase COLUMNS
}


void addElement(LockedArray_t* la, int row, step_t* element) {
	assert (row < la->rows);

	if (row == -1){
		row = rand() % la->rows;
	}

#ifdef PARALLEL
   pthread_mutex_lock(&la->locks[row]);
#endif

    int col = findEmptyColumn(la->arrays[row], la->columns);
	la->arrays[row][col] = element;  // Add the element

#ifdef PARALLEL
	pthread_mutex_unlock(&la->locks[row]);
#endif
}

void parallelDestroy(void* la_v, void* *helper, size_t length, size_t start, size_t end){
    LockedArray_t* la = (LockedArray_t*) la_v;
	for (int i = start; i < end; i++) {
		for (int j = 0; j < la->columns; j++){
			if (la->arrays[i][j] != NULL){
				step_free(la->arrays[i][j]);
			}
		}
        free(la->arrays[i]);
#ifdef PARALLEL
        pthread_mutex_destroy(&la->locks[i]);
#endif
    }
}

void destroyLockedArray(LockedArray_t* la) {
#ifdef PARALLEL
	parallel_for_c(la, NULL, 0, la->rows, BLOCKSIZE, parallelDestroy);
#else
	parallelDestroy(la, NULL, 0, 0, la->rows); // last argument was l
#endif
    free(la->arrays);
#ifdef PARALLEL
    free(la->locks);
#endif
    free(la);
}

/******************************************************************************/
/* END Lock Arrays                                                            */
/******************************************************************************/
#endif

/******************************************************************************/
/* Associative Smoother                                                       */
/******************************************************************************/


void buildFilteringElement(kalman_t* kalman, int i) {
	/*
		we denote by Z the matrix denoted by C in the article,
		because we already use C for the covariance of the
		observations.
	*/

	step_t* step_i = farray_get(kalman->steps, i);
	int n_i = step_i->dimension;
	int step = step_i->step;


	if (step == 0) {
		return;
	}

	if (step == 1) {
		step_t* step_1 = farray_get(kalman->steps, 0);
		matrix_t* G_i = step_1->G;
		matrix_t* o_i = step_1->o;
		matrix_t* C_i = step_1->C;

		char C_type = step_i->C_type;
						
		matrix_t* W_i_G_i = cov_weigh(C_i, C_type, G_i);
		matrix_t* W_i_o_i = cov_weigh(C_i, C_type, o_i);

		
		matrix_t* R = matrix_create_copy(W_i_G_i);
		matrix_t* Q = matrix_create_mutate_qr(R);
		
		matrix_mutate_apply_qt(R,Q,W_i_o_i);
		
		matrix_mutate_triu(R);
		matrix_t* m0 = matrix_create_trisolve(R,W_i_o_i);

		matrix_t* RT = matrix_create_transpose(R);
		matrix_t* RTR = matrix_create_multiply(RT,R);
		matrix_t* P0 = matrix_create_inverse(RTR);


		step_1->state = m0;
		step_1->covariance = P0;
		
		matrix_free(Q);
		matrix_free(R);

		matrix_free(W_i_G_i);
		matrix_free(W_i_o_i);
		
		matrix_free(RT);
		matrix_free(RTR);
	}

	matrix_t* F_i = step_i->F;
	matrix_t* c_i = step_i->c;
	matrix_t* K_i = explicit(step_i->K,step_i->K_type);

	if (step == 1) {
		step_t* step_1 = farray_get(kalman->steps, 0);
		matrix_t* F_iT = matrix_create_transpose(F_i);
		matrix_t* P0 = step_1->covariance;
		matrix_t* FP0 = matrix_create_multiply(F_i,P0);
		matrix_t* FPFT = matrix_create_multiply(FP0,F_iT);
		matrix_t* new_K_i = matrix_create_add(K_i,FPFT);
		matrix_free(K_i);
		K_i = new_K_i;
		matrix_free(F_iT);
		matrix_free(FP0);
		matrix_free(FPFT);
	}

	if (step_i->o == NULL) {
		step_i->Z = K_i;
		if (step == 1) {
			step_i->A = matrix_create_constant(n_i,n_i,0.0);
			step_t* step_1 = farray_get(kalman->steps, 0);
			matrix_t* m0 = step_1->state;
			matrix_t* b = matrix_create_add(m0,c_i);
			step_i->b = b;
		}else{
			step_i->A = matrix_create_copy(F_i);	
			step_i->b = matrix_create_copy(c_i);
		}
		step_i->e = NULL;
		step_i->J = NULL;
	}else{ // there are observations
		matrix_t* G_i = step_i->G;
		matrix_t* o_i = step_i->o;
		matrix_t* C_i = explicit(step_i->C,step_i->C_type);

		matrix_t* G_iT = matrix_create_transpose(G_i);
		matrix_t* KGT = matrix_create_multiply(K_i,G_iT);
		matrix_t* GKGT = matrix_create_multiply(G_i,KGT);
		matrix_t* S = matrix_create_add(GKGT,C_i);

		matrix_free(G_iT);
		matrix_free(KGT);	
		matrix_free(GKGT);
		matrix_free(C_i);

		matrix_t* ST = matrix_create_transpose(S);

		matrix_t* G_i_trans_inv_S_T = matrix_create_mldivide(ST,G_i); 
		matrix_t* G_i_trans_inv_S = matrix_create_transpose(G_i_trans_inv_S_T);
		
		matrix_free(ST);
		matrix_free(G_i_trans_inv_S_T);
		
		matrix_t* K = matrix_create_multiply(K_i, G_i_trans_inv_S);

		if (step == 1) {
			step_i->A = matrix_create_constant(n_i,n_i,0.0);
			step_t* step_1 = farray_get(kalman->steps, 0);
			matrix_t* m0 = step_1->state;
			matrix_t* F_im = matrix_create_multiply(F_i,m0);
			matrix_t* m1 = matrix_create_add(F_im,c_i);
			matrix_t* G_im = matrix_create_multiply(G_i,m1);
			matrix_t* o_G_im = matrix_create_subtract(o_i,G_im);
			matrix_t* K_o_G_im = matrix_create_multiply(K,o_G_im);
			matrix_t* b = matrix_create_add(m1,K_o_G_im);
			step_i->b = b;

			matrix_t *KS = matrix_create_multiply(K,S);
			matrix_t* KT = matrix_create_transpose(K);
			matrix_t* KSKT = matrix_create_multiply(KS,KT);
			matrix_t* Z = matrix_create_subtract(K_i,KSKT);
			step_i->Z = Z;

			matrix_free(F_im);
			matrix_free(m1);
			matrix_free(G_im);
			matrix_free(o_G_im);
			matrix_free(K_o_G_im);

			matrix_free(KS);
			matrix_free(KT);
			matrix_free(KSKT);
		}else{
			matrix_t* GF = matrix_create_multiply(G_i,F_i);
			matrix_t* KGF = matrix_create_multiply(K,GF);
			matrix_t* A = matrix_create_subtract(F_i,KGF);
			step_i->A = A;

			matrix_free(GF);
			matrix_free(KGF);

			matrix_t* G_ic = matrix_create_multiply(G_i,c_i);
			matrix_t* o_G_ic = matrix_create_subtract(o_i,G_ic);
			matrix_t* K_o_G_ic = matrix_create_multiply(K,o_G_ic);
			matrix_t* b = matrix_create_add(c_i,K_o_G_ic);
			step_i->b = b;

			matrix_free(G_ic);
			matrix_free(o_G_ic);
			matrix_free(K_o_G_ic);

			matrix_t* KG = matrix_create_multiply(K,G_i);
			matrix_t* KGK_i = matrix_create_multiply(KG,K_i);
			matrix_t* Z = matrix_create_subtract(K_i,KGK_i);
			step_i->Z = Z;

			matrix_free(KG);
			matrix_free(KGK_i);
		}
		matrix_free(K_i);

		matrix_t* G_ic = matrix_create_multiply(G_i,c_i);
		matrix_t* o_G_ic = matrix_create_subtract(o_i,G_ic);
		matrix_t* FT = matrix_create_transpose(F_i);
		matrix_t* FTG = matrix_create_multiply(FT,G_i_trans_inv_S);
		matrix_t* e = matrix_create_multiply(FTG,o_G_ic);
		step_i->e = e;

		matrix_free(G_ic);
		matrix_free(o_G_ic);
		matrix_free(FT);
		
		matrix_t* GF = matrix_create_multiply(G_i,F_i);
		matrix_t* J = matrix_create_multiply(FTG,GF);
		step_i->J = J;

		matrix_free(FTG);
		matrix_free(GF);

		matrix_free(G_i_trans_inv_S);
		matrix_free(K);
		matrix_free(S);
	}
}

void buildSmoothingElement(kalman_t* kalman, int i) {

	if (i == farray_size(kalman->steps) - 1) {
		step_t* step_i = farray_get(kalman->steps,i);
		int ni = step_i->dimension;
		step_i->E = matrix_create_constant(ni,ni,0.0);
		step_i->g = matrix_create_copy(step_i->state);
		step_i->L = matrix_create_copy(step_i->covariance);
	}else{
		step_t* step_i = farray_get(kalman->steps,i);
		matrix_t* x = step_i->state;
		matrix_t* P = explicit(step_i->covariance,'C');
		step_t* step_ip1 = farray_get(kalman->steps,i+1);
		matrix_t* F = step_ip1->F;
		matrix_t* Q = explicit(step_ip1->K,step_ip1->K_type);
		matrix_t* c = step_ip1->c;

		matrix_t* FT = matrix_create_transpose(F);
		matrix_t* PFT = matrix_create_multiply(P,FT);
		matrix_t* FPFT = matrix_create_multiply(F,PFT);
		matrix_t* FPFT_Q = matrix_create_add(FPFT,Q);

		matrix_t* PFT_T = matrix_create_transpose(PFT); 
		matrix_t* FPFT_Q_T = matrix_create_transpose(FPFT_Q);
		matrix_t* E_T = matrix_create_mldivide(FPFT_Q_T,PFT_T);
		matrix_t* E = matrix_create_transpose(E_T);

		step_i->E = E;

		matrix_free(FT);
		matrix_free(PFT);
		matrix_free(FPFT);
		matrix_free(FPFT_Q);
		matrix_free(PFT_T);
		matrix_free(FPFT_Q_T);
		matrix_free(E_T);


		matrix_t* Fx = matrix_create_multiply(F,x);
		matrix_t* Fx_c = matrix_create_add(Fx,c);
		matrix_t* E_Fx_c = matrix_create_multiply(E,Fx_c);
		matrix_t* g = matrix_create_subtract(x,E_Fx_c);
		step_i->g = g;

		matrix_free(Fx);
		matrix_free(Fx_c);
		matrix_free(E_Fx_c);
		matrix_t* EF = matrix_create_multiply(E,F);
		matrix_t* EFP = matrix_create_multiply(EF,P);
		matrix_t* L = matrix_create_subtract(P,EFP);
		step_i->L = L;

		matrix_free(EF);
		matrix_free(EFP);
		
		matrix_free(P);
		matrix_free(Q);
	}
}

void buildFilteringElements(void* kalman_v, void* *helper, size_t l, size_t start, size_t end){
  kalman_t* kalman = (kalman_t*) kalman_v;
	for (int j=start; j < end; j++) {
		buildFilteringElement(kalman,j);
	}
}
void buildSmoothingElements(void* kalman_v, void* *helper, size_t l, size_t start, size_t end){
  kalman_t* kalman = (kalman_t*) kalman_v;
	
	for (int j=start; j < end; j++) {
		buildSmoothingElement(kalman,j);
	}
}

void filtered_to_state(void* kalman_v, void* *filtered_v, size_t l, size_t start, size_t end){
  kalman_t* kalman = (kalman_t*) kalman_v;
  step_t* *filtered = (step_t**) filtered_v;

	int i = start;
	for (int j = start + 1; j < end + 1; j++) {
		step_t* step_j = farray_get(kalman->steps,j);
		step_t* filtered_i = filtered[i];
		step_j->state = matrix_create_copy(filtered_i->b);
		step_j->covariance = matrix_create_copy(filtered_i->Z);		
		if (i != 0){
			step_free(filtered_i);
		}
		i++;
	}
}
void smoothed_to_state(void* kalman_v, void* *smoothed_v, size_t l, size_t start, size_t end){
  kalman_t* kalman = (kalman_t*) kalman_v;
  step_t* *smoothed = (step_t**) smoothed_v;
	int i = 0;
	for (int j = start; j < end; j++) {
		i = l - 2 - j + 1;
		step_t* step_j = farray_get(kalman->steps,j);
		step_t* smoothed_i = smoothed[i];
		matrix_free(step_j->state);
		step_j->state = matrix_create_copy(smoothed_i->g);
		matrix_free(step_j->covariance);
		step_j->covariance = matrix_create_copy(smoothed_i->L);
		step_free(smoothed_i);
	}
}
void smooth(kalman_t* kalman) {
	int l = farray_size(kalman->steps);

	
#ifdef PARALLEL
	parallel_for_c(kalman, NULL, l, l, BLOCKSIZE, buildFilteringElements);
#else
	buildFilteringElements(kalman, NULL, l, 0, l);
#endif



#ifdef PARALLEL
	step_t** filtered = (step_t**)malloc((l - 1) * sizeof(step_t*));
	LockedArray_t* filtered_created_steps = initLockedArray(l);
	parallel_scan_c(kalman->steps->elements, (void**) filtered, filtered_created_steps, filteringAssociativeOperation , l - 1, 1);
	destroyLockedArray(filtered_created_steps);
#else
	step_t** filtered = cummulativeSumsSequential(kalman, filteringAssociativeOperation, 1, l - 1, 1);
#endif
	

#ifdef PARALLEL
	parallel_for_c(kalman, (void**) filtered, l, l - 1, BLOCKSIZE, filtered_to_state);
#else
	filtered_to_state(kalman, (void**) filtered, l, 0, l - 1);
#endif

	free(filtered);

#ifdef PARALLEL
	parallel_for_c(kalman, NULL, l, l, BLOCKSIZE, buildSmoothingElements);
#else
	buildSmoothingElements(kalman, NULL, l, 0, l);
#endif

#ifdef PARALLEL
	step_t** smoothed = (step_t**)malloc(l * sizeof(step_t*));
	LockedArray_t* smoothed_created_steps = initLockedArray(l);
	parallel_scan_c(kalman->steps->elements, (void**) smoothed, smoothed_created_steps, smoothingAssociativeOperation , l, -1);
	destroyLockedArray(smoothed_created_steps);
#else
	step_t** smoothed = cummulativeSumsSequential(kalman, smoothingAssociativeOperation, l - 1, 0, -1);
#endif

#ifdef PARALLEL
	parallel_for_c(kalman, (void**) smoothed, l, l - 1, BLOCKSIZE, smoothed_to_state);
#else
	smoothed_to_state(kalman, (void**) smoothed, l, 0, l - 1);
#endif

	free(smoothed);
}

//step_t** cummulativeSumsSequential(kalman_t* kalman, step_t* (*f)(step_t*, step_t*), int s, int e, int stride) {
step_t** cummulativeSumsSequential(kalman_t* kalman, void* (*f)(void*, void*, void*, int, int), int s, int e, int stride) {
		step_t** sums = (step_t**)malloc((abs(e-s) + 1)*sizeof(step_t*));
		int i = 0;
		step_t** a = (step_t**)kalman->steps->elements;
		step_t* sum = a[s];
		sums[i] = sum;
		i++;

		if (s > e) {
			for (int j=s+stride; j>=e; j+=stride) {
			        sum = (step_t*) f(sum,a[j], NULL, -1, -1);
				sums[i] = sum;
				i++;
			}
		}else{
			for (int j=s+stride; j<=e; j+=stride) {
				sum = (step_t*) f(sum,a[j], NULL, -1, -1);
				sums[i] = sum;
				i++;
			}
		}
		return sums;
}

//step_t* filteringAssociativeOperation(step_t* si, step_t* sj, LockedArray_t* created_steps, int row, int is_final_scan) {
void* filteringAssociativeOperation(void* si_v, void* sj_v, void* created_steps_v, int row, int is_final_scan) {
  step_t* si = (step_t*) si_v;
  step_t* sj = (step_t*) sj_v;
  LockedArray_t* created_steps = (LockedArray_t*) created_steps_v;
  
		if (si == NULL){
			return sj;
		}

		if (sj == NULL){
			return si;
		}


		step_t* sij = step_create();
#ifdef PARALLEL
		if (!is_final_scan){
			addElement(created_steps, row, sij);
		}
#endif
		int ni = matrix_rows(si->b);

		matrix_t* eye_ni = matrix_create_identity(ni,ni);
		matrix_t* siZ_sjJ = matrix_create_multiply(si->Z,sj->J);
		matrix_t* eye_ni_plus_siZ_sjJ = matrix_create_add(eye_ni,siZ_sjJ);

		matrix_t* AT = matrix_create_transpose(sj->A);
		matrix_t* other_T = matrix_create_transpose(eye_ni_plus_siZ_sjJ);
		matrix_t* XT = matrix_create_mldivide(other_T,AT);
		matrix_t* X = matrix_create_transpose(XT);
		
		
		matrix_free(siZ_sjJ);
		matrix_free(eye_ni_plus_siZ_sjJ);
		matrix_free(AT);
		matrix_free(other_T);
		matrix_free(XT);


		matrix_t* sjJ_siZ = matrix_create_multiply(sj->J,si->Z);
		matrix_t* eye_ni_plus_sjJ_siZ = matrix_create_add(eye_ni,sjJ_siZ);
		matrix_t* A_iT = matrix_create_transpose(si->A);

		
		other_T = matrix_create_transpose(eye_ni_plus_sjJ_siZ);
		matrix_t* YT = matrix_create_mldivide(other_T,si->A);
		matrix_t* Y = matrix_create_transpose(YT);

		matrix_free(sjJ_siZ);
		matrix_free(eye_ni_plus_sjJ_siZ);
		matrix_free(eye_ni);
		matrix_free(A_iT);
		matrix_free(other_T);
		matrix_free(YT);

		matrix_t* XA = matrix_create_multiply(X,si->A);
		sij->A = XA;

		matrix_t* siZ_sj_e = matrix_create_multiply(si->Z,sj->e);
		matrix_t* siZ_sj_e_plus_si_b = matrix_create_add(siZ_sj_e,si->b);
		matrix_t* X_siZ_sj_e_plus_si_b = matrix_create_multiply(X,siZ_sj_e_plus_si_b);
		matrix_t* b = matrix_create_add(X_siZ_sj_e_plus_si_b,sj->b);
		sij->b = b;

		matrix_free(siZ_sj_e);
		matrix_free(siZ_sj_e_plus_si_b);
		matrix_free(X_siZ_sj_e_plus_si_b);

		matrix_t* X_siZ = matrix_create_multiply(X,si->Z);
		matrix_t* A_jT = matrix_create_transpose(sj->A);
		matrix_t* X_siZ_AT = matrix_create_multiply(X_siZ,A_jT);
		matrix_t* Z = matrix_create_add(X_siZ_AT,sj->Z);
		sij->Z = Z;

		matrix_free(X_siZ);
		matrix_free(A_jT);
		matrix_free(X_siZ_AT);

		matrix_t* sjJ_si_b = matrix_create_multiply(sj->J,si->b);
		matrix_t* sj_e_minus_sjJ_si_b = matrix_create_subtract(sj->e,sjJ_si_b);
		matrix_t* Y_sj_e_minus_sjJ_si_b = matrix_create_multiply(Y,sj_e_minus_sjJ_si_b);
		matrix_t* e = matrix_create_add(Y_sj_e_minus_sjJ_si_b,si->e);
		sij->e = e;

		matrix_free(sjJ_si_b);
		matrix_free(sj_e_minus_sjJ_si_b);
		matrix_free(Y_sj_e_minus_sjJ_si_b);

		matrix_t* sjJ_siA = matrix_create_multiply(sj->J,si->A);
		matrix_t* Y_sjJ_siA = matrix_create_multiply(Y,sjJ_siA);
		matrix_t* J = matrix_create_add(Y_sjJ_siA,si->J);
		sij->J = J;

		matrix_free(sjJ_siA);
		matrix_free(Y_sjJ_siA);
		matrix_free(X);
		matrix_free(Y);

		return sij;

}


//step_t* smoothingAssociativeOperation(step_t* si, step_t* sj, LockedArray_t* created_steps, int row, int is_final_scan) {
void* smoothingAssociativeOperation(void* si_v, void* sj_v, void* created_steps_v, int row, int is_final_scan) {
  step_t* si = (step_t*) si_v;
  step_t* sj = (step_t*) sj_v;
  LockedArray_t* created_steps = (LockedArray_t*) created_steps_v;
		if (si == NULL){
			return sj;
		}

		if (sj == NULL){
			return si;
		}

		step_t* sij = step_create();
#ifdef PARALLEL
		if (!is_final_scan){
			addElement(created_steps, row, sij);
		}
#endif
		
		matrix_t* E = matrix_create_multiply(sj->E,si->E);
		sij->E = E;

		matrix_t* Eg = matrix_create_multiply(sj->E,si->g);
		matrix_t* g = matrix_create_add(Eg,sj->g);
		sij->g = g;

		matrix_free(Eg);

		matrix_t* ET = matrix_create_transpose(sj->E);
		matrix_t* EL = matrix_create_multiply(sj->E,si->L);
		matrix_t* ELT = matrix_create_multiply(EL,ET);
		matrix_t* L = matrix_create_add(ELT,sj->L);
		sij->L = L;

		matrix_free(ET);
		matrix_free(EL);
		matrix_free(ELT);

		return sij;
}

/******************************************************************************/
/* END Associative Smoother                                                   */
/******************************************************************************/

void apply_Q_on_block_matrix (matrix_t* R, matrix_t* Q, matrix_t** upper, matrix_t** lower) {
	matrix_t* concat = matrix_create_vconcat(*upper, *lower);
	matrix_mutate_apply_qt(R, Q, concat);

	matrix_t *new_upper = matrix_create_sub(concat,0,matrix_rows(*upper), 0, matrix_cols(concat));
	matrix_t *new_lower = matrix_create_sub(concat,matrix_rows(*upper),matrix_rows(*lower), 0, matrix_cols(concat));

	matrix_free(*upper);
	matrix_free(*lower);
	matrix_free(concat);

	*upper = new_upper;
	*lower = new_lower;
}
void assign_and_free(matrix_t ** original, matrix_t* new){
	if (*original != NULL){
		free(*original);
	}
	*original = new;
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

				matrix_free(step->G);
				matrix_free(step->o);
				matrix_free(step->C);

				matrix_free(step->Z);

				matrix_free(step->A);
				matrix_free(step->b);

				matrix_free(step->e);
				matrix_free(step->J);

				matrix_free(step->E);
				matrix_free(step->g);
				matrix_free(step->L);
				
				matrix_free(step->state);
				matrix_free(step->covariance);
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
/* MAIN                                                                       */
/******************************************************************************/


#if 0
int main(int argc, char *argv[]) {
	printf("hello world\n");
}


#if 0
  matrix_t* A = matrix_create(7,8);
  matrix_t* I = matrix_create_identity(7,8);
  matrix_set(I,6,7,13.0);
  if (debug) matrix_print(I,NULL,NULL);
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
  if (debug) matrix_print(W,NULL,NULL);
}
#endif

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
