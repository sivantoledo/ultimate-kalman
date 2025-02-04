/*
 * ultimatekalman.c
 *
 * (C) Sivan Toledo, 2022
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// for getenv
#include <stdlib.h>

//#include "parallel_for_c.h"
void parallel_for_c(void* kalman, void* indices, int length, int** helper, size_t n, size_t block_size, void (*func)(void*, void*, int, int**, size_t, size_t));
/*
void assign_indices(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end) {
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**) indices_v;
	for (int i = start; i < end; ++i) {
		indices[i] = farray_get(kalman->steps,i);
	}
}

void f(kalman_t* kalman) {

	int length = farray_size(kalman->steps);
	step_t* *indices = (step_t**)malloc(length * sizeof(step_t*));

	parallel_for_c(kalman, indices, length, NULL, length, BLOCKSIZE, assign_indices);
}
*/
int kalman_parallel_init(int number_of_threads);
int kalman_parallel_blocksize(int blocksize_in);

//#define BLOCKSIZE 1000
#define BLOCKSIZE 10

#if defined(PARALLEL) && !defined(MACOS)
#define malloc(x) aligned_alloc(64,(x))
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
#include "ultimatekalman_oddeven.h"

char kalman_covariance_type(kalman_t* kalman, int64_t si) { return 'C'; }


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
	//fprintf(stderr,"sizeof matrix_t is %lu\n",sizeof(matrix_t));
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

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

static step_t* step_create() {
	step_t* s = malloc(sizeof(step_t));
	s->step      = -1;
	s->dimension = -1;


	s->H = NULL;
	s->F = NULL;
	s->c = NULL;

	s->G = NULL;
	s->o = NULL;
	s->C = NULL;

	s->X = NULL;
	s->Y = NULL;
	s->Z = NULL;
	s->R = NULL;

	s->X_tilde = NULL;
	s->F_tilde = NULL;
	s->H_tilde = NULL;
	s->R_tilde = NULL;
	s->G_tilde = NULL;

	s->state = NULL;
	s->covariance = NULL;

	assert( s != NULL );
	return s;
}

void step_free(step_t* s) {
	matrix_free(s->H);
	matrix_free(s->F);
	matrix_free(s->c);

	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->X);
	matrix_free(s->Y);
	matrix_free(s->Z);
	matrix_free(s->R);

	matrix_free(s->X_tilde);
	matrix_free(s->F_tilde);
	matrix_free(s->H_tilde);
	matrix_free(s->R_tilde);
	matrix_free(s->G_tilde);

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

	char* blocksize_string = getenv("TBB_BLOCKSIZE");
	int blocksize = 0;
	if (blocksize_string != NULL && sscanf(blocksize_string,"%d",&blocksize)==1) {
		kalman_parallel_blocksize(blocksize);
		printf("setting blocksize to %d\n",blocksize);
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


	matrix_t* V_i_H_i = cov_weigh(K_i,K_type,H_i);
	matrix_t* V_i_F_i = cov_weigh(K_i,K_type,F_i);
	matrix_t* V_i_c_i = cov_weigh(K_i,K_type,c_i);

	matrix_mutate_scale(V_i_F_i,-1.0);

	kalman->current->H	= V_i_H_i;
	kalman->current->F 	= V_i_F_i;
	kalman->current->c  = V_i_c_i;
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

	matrix_t* W_i_G_i = NULL;
	matrix_t* W_i_o_i = NULL;

	if (o_i != NULL) {
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("C_i(%c) ",C_type);
		matrix_print(C_i,"%.3e");
		printf("G_i ");
		matrix_print(G_i,"%.3e");
		printf("o_i ");
		matrix_print(o_i,"%.3e");
#endif
	W_i_G_i = cov_weigh(C_i,C_type,G_i);
	W_i_o_i = cov_weigh(C_i,C_type,o_i);

	kalman->current->G  = W_i_G_i;
	kalman->current->o  = W_i_o_i;

#ifdef BUILD_DEBUG_PRINTOUTS
		printf("W_i_G_i ");
		matrix_print(W_i_G_i,"%.3e");
		printf("W_i_o_i ");
		matrix_print(W_i_o_i,"%.3e");
#endif
	}

	farray_append(kalman->steps,kalman->current);
	kalman->current = NULL;

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

void assign_indices(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end) {
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**) indices_v;
	for (int i = start; i < end; ++i) {
		indices[i] = farray_get(kalman->steps,i);
	}
}

void kalman_smooth(kalman_t* kalman) {

	int length = farray_size(kalman->steps);
	step_t* *indices = (step_t**)malloc(length * sizeof(step_t*));

#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, length, BLOCKSIZE, assign_indices);
#else
	assign_indices(kalman, indices, length, NULL, 0, length);
#endif

	//for (int i = 0; i < length; ++i) {
	//	indices[i] = farray_get(kalman->steps,i);
	//}

	parallel_smooth(kalman, indices, length);

	free(indices);
}

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
void free_and_assign(matrix_t** original, matrix_t* new){
	if (*original != NULL){
		matrix_free(*original);
	}
	*original = new;
}
void G_F_to_R_tilde(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; ++j_) {
		int j = j_ * 2;

		step_t* step_i = indices[j];
		matrix_t* G_i = step_i->G;

		if (j == length - 1){
			matrix_t* o_i = step_i->o;

			matrix_t* R_tilde = matrix_create_copy(G_i);

			matrix_t* TAU = matrix_create_mutate_qr(R_tilde);
			matrix_mutate_apply_qt(R_tilde,TAU,o_i);

			step_i->R_tilde = R_tilde;
			matrix_mutate_triu(step_i->R_tilde);
			matrix_free(TAU);
		}else{

			step_t* step_ipo = indices[j + 1];
			matrix_t* o_i = step_i->o;

			matrix_t* c_ipo = step_ipo->c;

			matrix_t* F_ipo = step_ipo->F;
			matrix_t* H_ipo = step_ipo->H;

			matrix_t* concat = matrix_create_vconcat(G_i, F_ipo);

			matrix_t* Q, *R_tilde;

			R_tilde = concat;
			Q = matrix_create_mutate_qr(R_tilde);

			step_i->R_tilde = matrix_create_sub(R_tilde,0,matrix_cols(R_tilde), 0, matrix_cols(R_tilde));
			matrix_mutate_triu(step_i->R_tilde);
			step_i->X = matrix_create_constant(matrix_rows(G_i), matrix_cols(H_ipo), 0);

			free_and_assign(&(step_ipo->H_tilde), matrix_create_copy(H_ipo));

			apply_Q_on_block_matrix(R_tilde, Q, &(step_i->o), &(step_ipo->c));
			apply_Q_on_block_matrix(R_tilde, Q, &(step_i->X), &(step_ipo->H_tilde));
			matrix_free(R_tilde);
			matrix_free(Q);
		}
	}//Done first part of the algorithm
}
void H_R_tilde_to_R(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; j_++){
		int j = j_ * 2;

		step_t* step_i = indices[j];
		if (j == 0){ //First index
			step_i->R = matrix_create_copy(step_i->R_tilde);
			continue;
		}

		matrix_t* R_tilde = step_i->R_tilde;
		matrix_t* F_i = step_i->F;
		matrix_t* H_i = step_i->H;

		matrix_t* concat = matrix_create_vconcat(H_i, R_tilde);
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));
		
		
		matrix_mutate_triu(R);
		step_i->R = R;

		step_i->Z = matrix_create_constant(matrix_rows(R_tilde), matrix_cols(F_i), 0);

		step_i->F_tilde = matrix_create_copy(F_i);
		apply_Q_on_block_matrix(R_tall, Q, &(step_i->F_tilde), &(step_i->Z));
		apply_Q_on_block_matrix(R_tall, Q, &(step_i->c), &(step_i->o));

		if(j + 1 != length){
			matrix_t *X = step_i->X;

			step_i->Y = matrix_create_constant(matrix_rows(F_i), matrix_cols(X), 0);

			step_i->X_tilde = matrix_create_copy(X);

			apply_Q_on_block_matrix(R_tall, Q, &(step_i->Y), &(step_i->X_tilde));

		}
		matrix_free(R_tall);
		matrix_free(Q);
	} //Done second part of the algorithm

}
void H_tilde_G_to_G_tilde(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; j_++){
		int j = j_ * 2;
		// Sivan Feb 2025 wrong type and seens not to be used, commenting out
		//int i = indices[j];

		step_t* step_ipo = indices[j + 1];
		matrix_t* H_tilde = step_ipo->H_tilde;
		matrix_t* G_ipo = step_ipo->G;

		matrix_t* concat = matrix_create_vconcat(H_tilde, G_ipo);
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));

		matrix_mutate_triu(R);

		free_and_assign(&(step_ipo->G_tilde), R); 


		apply_Q_on_block_matrix(R_tall, Q, &(step_ipo->c), &(step_ipo->o));
		matrix_free(R_tall);
		matrix_free(Q);


	}//Done third part of the algorithm
}
void Variables_Renaming(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; ++j_){
		int j = j_ * 2;
		// int i = indices[j];
		// int ipo = indices[j + 1];
		step_t* step_ipo = indices[j + 1];
		matrix_t* G_tilde = step_ipo->G_tilde;
		matrix_t* o = step_ipo->c;

		matrix_t * copy_G_tilde = matrix_create_copy(G_tilde);
		matrix_t * copy_o = matrix_create_copy(o);

		if (j != 0){
			step_t* step_i = indices[j];
			matrix_t* Z = step_i->Z;
			matrix_t* X_tilde = step_i->X_tilde;
			matrix_t* c = step_i->o;

			matrix_t * copy_Z = matrix_create_copy(Z);
			matrix_t * copy_X_tilde = matrix_create_copy(X_tilde);
			matrix_t * copy_c = matrix_create_copy(c);

			free_and_assign(&(step_ipo->F), copy_Z); 

			free_and_assign(&(step_ipo->H), copy_X_tilde); 

			free_and_assign(&(step_ipo->c), copy_c);
		}
		free_and_assign(&(step_ipo->G), copy_G_tilde);
		free_and_assign(&(step_ipo->o), copy_o);
	}
}

void Init_new_indices(void* new_indices_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	step_t** new_indices = (step_t**) new_indices_v;
	step_t** indices     = (step_t**) indices_v;

    for (int i = start; i < end; ++i) {
		new_indices[i] = indices[2*i + 1];
	}
}
void Solve_Estimates(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; ++j_){
		int j = j_ * 2;

		step_t* step_i = indices[j];
		matrix_t* R = step_i->R;

		if (j == 0){

			step_t* step_ipo = indices[j + 1];
			matrix_t* x_ipo = step_ipo->state;

			matrix_t* X = step_i->X;
			matrix_t* o_i = step_i->o;
			matrix_t* mul = matrix_create_constant(matrix_rows(X), 1, 0);
			matrix_mutate_gemm(1, X, x_ipo, 0, mul);
			matrix_t* new_b = matrix_mutate_subtract(o_i, mul);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul);
			matrix_free(new_b);

		}else if (j != length - 1){
			
			step_t* step_ipo = indices[j + 1];
			matrix_t* x_ipo = step_ipo->state;
			step_t* step_imo = indices[j - 1];
			matrix_t* x_imo = step_imo->state;

			matrix_t* F_tilde = step_i->F_tilde;
			matrix_t* Y = step_i->Y;
			matrix_t* c = step_i->c;

			matrix_t* mul1 = matrix_create_constant(matrix_rows(F_tilde), 1, 0);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul1);
			matrix_t* mul2 = matrix_create_constant(matrix_rows(Y), 1, 0);
			matrix_mutate_gemm(1, Y, x_ipo, 0, mul2);
			matrix_t* new_b_mid = matrix_mutate_subtract(c, mul1);
			matrix_t* new_b = matrix_mutate_subtract(new_b_mid, mul2);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul1);
			matrix_free(mul2);
			matrix_free(new_b_mid);
			matrix_free(new_b);
		}else{

			
			step_t* step_imo = indices[j - 1];
			matrix_t* x_imo = step_imo->state;

			matrix_t* F_tilde = step_i->F_tilde;
			matrix_t* c = step_i->c;

			matrix_t* mul = matrix_create_constant(matrix_rows(F_tilde), 1, 0);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul);
			matrix_t* new_b = matrix_mutate_subtract(c, mul);

			matrix_mutate_triu(R);

			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul);
			matrix_free(new_b);
		}
	}
}
// ==========================================
// Cov Change 4
// ==========================================

void Convert_LDLT(void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	for (int j_ = start; j_ < end; ++j_){
		int j = j_ * 2;

		step_t* step = indices[j];
		matrix_t* R = step->R;
		matrix_t* R_inv = matrix_create_inverse(R);
		step->R = R_inv;
		
		if (j == 0){
			matrix_t* X = step->X;
			matrix_t* X_inv = matrix_create_trisolve(R,X);
			matrix_free(X);
			step->X = X_inv;

		}else if (j != length - 1){
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve(R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;	

			matrix_t* Y = step->Y;
			matrix_t* Y_inv = matrix_create_trisolve(R,Y);
			matrix_free(Y);
			step->Y = Y_inv;	
			
		}else{
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve(R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;		
		}
	}
}

void SelInv(void* kalman_v, void* indices_v, int length, int** converters, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;

	int * vtop_n = converters[0];
	int * ptov_n2 = converters[1];
	
	for (int j_ = start; j_ < end; ++j_){
		int j = j_ * 2;

		step_t* step = indices[j];
	// 	Ainv(inz,j) = -Ainv(inz,inz) * L(inz,j);
	// 	Ainv(j,inz) = Ainv(inz,j)';
	// 	Ainv(j,j) = 1/D(j,j) - Ainv(j,inz) * L(inz,j);
		
		if (j == 0){
			int inz1 = 0;
			
			int phisical_column = ceil(((double)length)/2) + ptov_n2[inz1];
			matrix_t* L_j_inz = step->X;

			step_t *phisical_step = indices[vtop_n[phisical_column]];
			int phisical_row = phisical_step->step;

			matrix_t* Ainv_inz_inz = phisical_step->R;

			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
			matrix_t* Dinv = step->R;

			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_mutate_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);
			matrix_free(Dinv_Dinv_t);

			matrix_free(step->X);
			step->X = Ainv_j_inz;

			matrix_free(step->R);
			step->R = Ainv_j_j;

		}else if (j != length - 1){
			int inz1 = j/2 - 1;
			int inz2 = j/2;

			int phisical_column1 = ceil(((double)length)/2) + ptov_n2[inz1];

			int phisical_column2 = ceil(((double)length)/2) + ptov_n2[inz2];

			int regular_order = 1;

	// 		We have two cases here: Y,F_tilde and F_tilde,Y 
	// 		and the order is matter.

			step_t *phisical_step1;
			step_t *phisical_step2;
			
			matrix_t* L_j_inz;

			if (phisical_column1 < phisical_column2){ //This is the ragular order

				matrix_t* F_tilde_t = matrix_create_transpose(step->F_tilde);
				matrix_t* Y_t = matrix_create_transpose(step->Y);
				matrix_t* L_j_inz_t = matrix_create_vconcat(F_tilde_t, Y_t);
				L_j_inz = matrix_create_transpose(L_j_inz_t);

				matrix_free(F_tilde_t);
				matrix_free(Y_t);
				matrix_free(L_j_inz_t);

				phisical_step1 = indices[vtop_n[phisical_column1]];
				int phisical_row1 = phisical_step1->step;
				phisical_step2 = indices[vtop_n[phisical_column2]];
				int phisical_row2 = phisical_step2->step;
			}else{
				regular_order = 0;

				matrix_t* F_tilde_t = matrix_create_transpose(step->F_tilde);
				matrix_t* Y_t = matrix_create_transpose(step->Y);
				matrix_t* L_j_inz_t = matrix_create_vconcat(Y_t, F_tilde_t);
				L_j_inz = matrix_create_transpose(L_j_inz_t);

				matrix_free(F_tilde_t);
				matrix_free(Y_t);
				matrix_free(L_j_inz_t);
					
				int tmp = phisical_column1;
				phisical_column1 = phisical_column2;
				phisical_column2 = tmp;

				phisical_step1 = indices[vtop_n[phisical_column1]];
				int phisical_row1 = phisical_step1->step;
				phisical_step2 = indices[vtop_n[phisical_column2]];
				int phisical_row2 = phisical_step2->step;	
			}

			
	// 		Now we have to do the same process in row 2 but
	// 		column 2. but we can't call virtual_phisical

			int start_index = ceil(((double)length)/2);

			matrix_t* top_right;
	// 		again we have 3 cases:
			if (start_index == phisical_column1){
				top_right = phisical_step1->X;
			}else if (start_index + ceil(floor(((double)length)/2)/2) > phisical_column1 + 1 || (length/2)%2 == 0){
	
				int source_index = floor(((double)vtop_n[phisical_column2])/4);

	// 			Now lets check the source of F_tilde and Y of this row
				int F_tilde_location = phisical_column1 - start_index - 1;
				if (source_index == F_tilde_location){	
					top_right = phisical_step1->F_tilde;
				}else{
					top_right = phisical_step1->Y;
				    assert(source_index == phisical_column1 - start_index);
				}
			}else{	
				top_right = phisical_step1->F_tilde;
			}



	// 		Ainv_inz_inz = [kalman.steps{phisical_row1 + 1}.R , top_right
	// 						top_right'                        , kalman.steps{phisical_row2 + 1}.R];
			matrix_t* top_right_t = matrix_create_transpose(top_right);
			matrix_t* R_top_right_t = matrix_create_vconcat(phisical_step1->R, top_right_t);
			matrix_t* top_right_R = matrix_create_vconcat(top_right, phisical_step2->R);
			matrix_t* t1 = matrix_create_transpose(R_top_right_t);
			matrix_t* t2 = matrix_create_transpose(top_right_R);
			matrix_t* Ainv_inz_inz_t = matrix_create_vconcat(t1, t2);
			matrix_t* Ainv_inz_inz = matrix_create_transpose(Ainv_inz_inz_t);

			matrix_free(top_right_t);
			matrix_free(R_top_right_t);
			matrix_free(top_right_R);
			matrix_free(t1);
			matrix_free(t2);
			matrix_free(Ainv_inz_inz_t);

	// 		Ainv_j_inz = -L_j_inz * Ainv_inz_inz;
			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_mutate_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);
			
			if (regular_order){ //This is the ragular order

				int F_tilde_cols = matrix_cols(step->F_tilde);
				int Y_cols = matrix_cols(step->Y);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, F_tilde_cols);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), F_tilde_cols, matrix_cols(Ainv_j_inz) - F_tilde_cols);	
			}else{
				int F_tilde_cols = matrix_cols(step->F_tilde);
				int Y_cols = matrix_cols(step->Y);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, Y_cols);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), Y_cols, matrix_cols(Ainv_j_inz) - Y_cols);
			}
			
			matrix_free(step->R);
			step->R = Ainv_j_j;

		}else{
			int inz1 = j/2 - 1;			

			int phisical_column = ceil(((double)length)/2) + ptov_n2[inz1];

			matrix_t* L_j_inz = step->F_tilde;

			step_t *phisical_step = indices[vtop_n[phisical_column]];
			int phisical_row = phisical_step->step;

			matrix_t* Ainv_inz_inz = phisical_step->R;

	// 		Ainv_j_inz = -L_j_inz * Ainv_inz_inz;
			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_mutate_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);

	// 		kalman.steps{i}.F_tilde = Ainv_j_inz;
			matrix_free(step->F_tilde);
			step->F_tilde = Ainv_j_inz;
	// 		kalman.steps{i}.R = Ainv_j_j;
			matrix_free(step->R);
			step->R = Ainv_j_j;
		}
		step->covariance = matrix_create_copy(step->R);
	}
}

// ==========================================
// Cov End Change 4
// ==========================================

// ==========================================
// Cov Change 2
// ==========================================


void one_layer_converter (void* kalman_v, void* indices_v, int length, int** helper, size_t start, size_t end){
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *indices = (step_t**)indices_v;
	 
	int* virtual_to_phisical = helper[0];
	int* phisical_to_virtual = helper[1];
	int* variables = helper[2];

	int n = variables[0];
	int start_ = variables[1];
	int jump = variables[2];
	int skip = variables[3];
	int index = variables[4];

	for (int j_ = start; j_ < end; ++j_){
		int i = j_ * jump + start_;

		virtual_to_phisical[index + j_] = i;
		if (!skip){
			phisical_to_virtual[i/2] = (index + j_) - ceil(((double)n)/2);
		}
	}
	variables[4] = index;
}

int ** virtual_phisical(int n){

	int* virtual_to_phisical = (int *)malloc(n * sizeof(int));
	int* phisical_to_virtual = (int *)malloc((n/2) * sizeof(int));
	
	int* variables = (int *)malloc(5 * sizeof(int));

	int start = 0;
	int jump = 2;

	int index = 0;

	int skip = 1;

	int ** result = (int **)malloc(3 * sizeof(int *));
	result[0] = virtual_to_phisical;
	result[1] = phisical_to_virtual;
	result[2] = variables;

	while (start < n){

		variables[0] = n;
		variables[1] = start;
		variables[2] = jump;
		variables[3] = skip;
		variables[4] = index;

		#ifdef PARALLEL
		parallel_for_c(NULL, NULL, 0, result, ceil((double)(n - start)/jump),BLOCKSIZE,one_layer_converter);
		#else
		one_layer_converter(NULL, NULL, 0, result, 0, ceil((double)(n - start)/jump));
		#endif

		index += ceil((double)(n - start)/jump);

		skip = 0;

		start = start * 2 + 1;
		jump = jump * 2;
	}
	
	free(variables);
	return result;
}
// ==========================================
// End Cov Change 2
// ==========================================

void parallel_smooth(kalman_t* kalman, step_t* *indices, int length) {


	
    if (length == 1) {

		step_t* singleStep = indices[0];

		matrix_t* G = matrix_create_copy(singleStep->G);
		matrix_t* o = matrix_create_copy(singleStep->o);

		matrix_t* TAU = matrix_create_mutate_qr(G);
		matrix_mutate_apply_qt(G,TAU,o);

		matrix_mutate_triu(G);

		// ==========================================
		// Cov Change 3
		// ==========================================

		matrix_t* RT = matrix_create_transpose(G);

		matrix_t* RT_R = matrix_create_constant(matrix_cols(G), matrix_cols(G), 0);
		matrix_mutate_gemm(1, RT, G, 0, RT_R);

		singleStep->R = matrix_create_inverse(RT_R);

		singleStep->covariance = matrix_create_copy(singleStep->R);

		matrix_free(RT);
		matrix_free(RT_R);
		
		// ==========================================
		// End Cov Change 3
		// ==========================================
		singleStep->state = matrix_create_trisolve(G,o);
		
		matrix_free(TAU);
		matrix_free(G);
		matrix_free(o);
		
		return;
    }
	// First part of the algorithm
	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, (length + 1)/2,BLOCKSIZE,G_F_to_R_tilde);
	#else
	G_F_to_R_tilde(kalman, indices, length, NULL, 0, (length + 1)/2);
	#endif

	//Second part of the algorithm
	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, (length + 1)/2,BLOCKSIZE, H_R_tilde_to_R);
	#else
	H_R_tilde_to_R(kalman, indices, length, NULL, 0, (length + 1)/2);
	#endif

	//Third part of the algorithm
	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, length/2,BLOCKSIZE, H_tilde_G_to_G_tilde);
	#else
	H_tilde_G_to_G_tilde(kalman, indices, length, NULL, 0, length/2);
	#endif

	//Fix the last index

	if (length % 2 == 1){
		step_t* step_lmo = indices[length - 2];
		matrix_t* G_tilde = step_lmo->G_tilde;
		matrix_t* o_i = step_lmo->c;

		step_t* step_last = indices[length - 1];
		matrix_t* Z = step_last->Z;
		
		matrix_t* c_i = step_last->o;
		matrix_t* concat = matrix_create_vconcat(G_tilde, Z);
		
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));

		matrix_mutate_triu(R);
		free_and_assign(&(step_lmo->G_tilde), R); //step_lmo->G_tilde = R;

		apply_Q_on_block_matrix(R_tall, Q, &(step_lmo->c), &(step_last->o));
		matrix_free(R_tall);
		matrix_free(Q);
	}

	// Create new Kalman for the recursion
	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, length/2, BLOCKSIZE, Variables_Renaming);
	#else
	Variables_Renaming(kalman, indices, length, NULL, 0, length/2);
	#endif
	

	// The Recursion

	step_t* *new_indices = (step_t**)malloc(length/2 * sizeof(step_t*));

	#ifdef PARALLEL
	parallel_for_c(new_indices, indices, length, NULL, length/2,BLOCKSIZE, Init_new_indices);
	#else
	Init_new_indices(new_indices, indices, length, NULL, 0, length/2);
	#endif
	
	
	parallel_smooth(kalman, new_indices, length/2);

	free(new_indices);

	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, (length + 1)/2, BLOCKSIZE, Solve_Estimates);
	#else
	Solve_Estimates(kalman, indices, length, NULL, 0, (length + 1)/2);
	#endif

	// ==========================================
	// Change 1
	// ==========================================

	// Since the selinv algorithm works with LDL^T matrices, we
	// start with converting our RR^T to LDL^T

	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, NULL, (length + 1)/2, BLOCKSIZE, Convert_LDLT);
	#else
	Convert_LDLT(kalman, indices, length, NULL, 0, (length + 1)/2);
	#endif

	// Now we can start the selinv algrithm for out case
	
	int ** result = virtual_phisical(length);
	

	#ifdef PARALLEL
	parallel_for_c(kalman, indices, length, result, (length + 1)/2, BLOCKSIZE, SelInv);
	#else
	SelInv(kalman, indices, length, result, 0, (length + 1)/2);
	#endif
	
	free(result[0]);
	free(result[1]);
	free(result);
	// % ==========================================
	// % End Change 1
	// % ==========================================
	

}
	// //solve
	// matrix_mutate_triu(kalman->current->Rdiag);
	// state = matrix_create_trisolve(kalman->current->Rdiag,y);

	// //also importent
	// matrix_mutate_chop(B,MIN(matrix_rows(B),n_imo),matrix_cols(B));
	// //importent COMBO
	// matrix_t* TAU = matrix_create_mutate_qr(A);
	// matrix_mutate_apply_qt(A,TAU,B);
	// matrix_mutate_apply_qt(A,TAU,y);
	// matrix_free(TAU);



matrix_t* kalman_covariance(kalman_t* kalman, int64_t si) {
	//int32_t i,j;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_covariance\n");
#endif

	if (farray_size(kalman->steps) == 0) return NULL;

	if (si < 0) si = farray_last_index(kalman->steps);
	step_t* step = farray_get(kalman->steps,si);

	int32_t n_i = step->dimension;

	matrix_t* cov = NULL;

	if (step->covariance != NULL) {
		cov = matrix_create_copy(step->covariance);

		//printf("cov ");
		//matrix_print(cov,"%.3e");

		//for (j=0; j<n_i; j++) {
		//	for (i=j+1; i<n_i; i++) {
		//		matrix_set(cov,i,j,0.0);
		//	}
		//}
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
	//printf("rollback %d\n",si);
	if (farray_size(kalman->steps) == 0) return;

	if (si > farray_last_index (kalman->steps)) return; // we can roll  back even the last step (its observation)
	if (si < farray_first_index(kalman->steps)) return;

	step_t* step;
	do {
		step = farray_drop_last(kalman->steps);
	 	if (step->step == si) {
	 		//printf("rollback got to %d\n",si);

            matrix_free(step->G);
            matrix_free(step->o);
            matrix_free(step->C);

            matrix_free(step->X);
            matrix_free(step->Y);
            matrix_free(step->Z);
            matrix_free(step->R);

            matrix_free(step->X_tilde);
            matrix_free(step->F_tilde);
            matrix_free(step->H_tilde);
            matrix_free(step->R_tilde);
            matrix_free(step->G_tilde);

            matrix_free(step->state);
            matrix_free(step->covariance);

	 		kalman->current = step;

	 	} else {
	// 		//printf("rollback dropped step %d\n",step->step);
	 	}
	} while (step->step > si);
	// //printf("rollback to %d new latest %d\n",si,kalman_latest(kalman));
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
