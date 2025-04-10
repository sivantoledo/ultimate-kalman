#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef BUILD_MEX
#include "mex.h"
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
#undef blas_int_t
#define blas_int_t BUILD_BLAS_INT
#define HAS_BLAS_INT
#endif

#ifdef BUILD_MKL
#define HAS_BLAS_H
#define HAS_LAPACK_H
#undef blas_int_t
#define blas_int_t MKL_INT
#define HAS_BLAS_INT
#ifdef BUILD_BLAS_UNDERSCORE
#undef BUILD_BLAS_UNDERSCORE
#endif
#include <mkl.h>
#endif

#ifndef HAS_BLAS_INT
#undef blas_int_t
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

/******************************************************************************/
/* MATRICES                                                                   */
/******************************************************************************/

typedef struct matrix_st {
	int32_t row_dim;
	int32_t col_dim;
	int32_t ld;      // leading dimension
	double* elements;
}
kalman_matrix_t
#ifdef KALMAN_MATRIX_SHORT_TYPE
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

/*
 * Creates a zero matrix
 */
kalman_matrix_t* matrix_create(int32_t rows, int32_t cols);
kalman_matrix_t* matrix_create_constant(int32_t rows, int32_t cols, double c);
kalman_matrix_t* matrix_create_copy(kalman_matrix_t* A);
kalman_matrix_t* matrix_create_sub(kalman_matrix_t* A, int32_t first_row, int32_t rows, int32_t first_col, int32_t cols);
void             matrix_mutate_copy_sub(kalman_matrix_t* C, int32_t first_row, int32_t first_col, kalman_matrix_t* A);
void             matrix_free(kalman_matrix_t* A);

int32_t matrix_rows(kalman_matrix_t* A);
int32_t matrix_cols(kalman_matrix_t* A);
int32_t matrix_ld  (kalman_matrix_t* A);

/*
 * chops the matrix to the given number of rows and columns, from
 * the (0,0) element.
 */
//void matrix_mutate_chop(kalman_matrix_t* A, int32_t rows, int32_t cols);

/*
 * Creates an identity matrix (can be rectangular; main diagonal is 1, rest 0).
 */
kalman_matrix_t* matrix_create_identity(int32_t rows, int32_t cols);

/*
 * Creates a matrix from a rowwise C matrix.
 */
kalman_matrix_t* matrix_create_from_rowwise(double* naked, int32_t rows, int32_t cols);

#ifdef BUILD_MEX
kalman_matrix_t* matrix_create_from_mxarray(const mxArray* mx);
mxArray*  matrix_copy_to_mxarray    (kalman_matrix_t* A);
#endif


void matrix_print(kalman_matrix_t* A, char* format);

/*
 * Intended mostly for testing that the BLAS library is working and linked correctly
 */
void matrix_mutate_gemm                 (double ALPHA, kalman_matrix_t* A, kalman_matrix_t* B, double BETA, kalman_matrix_t* C);

void matrix_mutate_scale                (kalman_matrix_t* A, double s);
void matrix_mutate_triu                 (kalman_matrix_t* A);
void matrix_mutate_chop                 (kalman_matrix_t* A, int32_t rows, int32_t cols);
void matrix_mutate_copy                 (kalman_matrix_t* C, kalman_matrix_t* A);
void matrix_mutate_trisolve             (char* triangle, kalman_matrix_t* U, kalman_matrix_t* b);
void matrix_mutate_chol                 (kalman_matrix_t* L);

/*
 * LAPACK-type QR factorization.
 *
 * Mutates A so that it contains both R (upper triangle) and
 * the reflectors whose product is Q. Returns TAU, the scalar
 * factors of the reflectors (a vector). See documentation of
 * DGEQRF for more details.
 */
kalman_matrix_t* matrix_create_mutate_qr(kalman_matrix_t* A);
/*
 * Apply Q^T, computed by matrix_create_mutate_qr, to a matrix C.
 * The first argument QR contains both R and most of the representation
 * of Q. TAU is another part of the representation. See documentation of
 * DORMQR for more details.
 */
void matrix_mutate_apply_qt             (kalman_matrix_t* QR, kalman_matrix_t* TAU, kalman_matrix_t* C);


kalman_matrix_t* matrix_create_trisolve (char* triangle, kalman_matrix_t* U, kalman_matrix_t* b);
kalman_matrix_t* matrix_create_chol     (kalman_matrix_t* L);
kalman_matrix_t* matrix_create_vconcat  (kalman_matrix_t* A, kalman_matrix_t* B);

kalman_matrix_t* matrix_create_inverse  (kalman_matrix_t* A);
kalman_matrix_t* matrix_create_transpose(kalman_matrix_t* A);
kalman_matrix_t* matrix_create_copy     (kalman_matrix_t* A);

kalman_matrix_t* matrix_create_mldivide (kalman_matrix_t* A, kalman_matrix_t* B);
kalman_matrix_t* matrix_create_multiply (kalman_matrix_t* A, kalman_matrix_t* B);
kalman_matrix_t* matrix_create_subtract (kalman_matrix_t* A, kalman_matrix_t* B);
kalman_matrix_t* matrix_create_add      (kalman_matrix_t* A, kalman_matrix_t* B);


//kalman_matrix_t* matrix_create_vconcat(kalman_matrix_t* A, kalman_matrix_t* B);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif /* ifndef MATRIX_OPS_H */
