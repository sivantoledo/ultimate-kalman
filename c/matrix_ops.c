/*
 * matrix_ops.c
 *
 * (C) Sivan Toledo, 2022-2025
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
// for "unused" attribute
#define __attribute__(x)
#include <float.h>
// string.h for memcpy
#else
#include <unistd.h>
#endif

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

#define KALMAN_MATRIX_SHORT_TYPE
#include "matrix_ops.h"
#include "memory.h"

/******************************************************************************/
/* UTILITIES                                                                  */
/******************************************************************************/

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

// static double NaN = 0.0 / 0.0;

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

/* function from ultimatekalman_oddeven.c; name is wrong, there is no mutation */
/*
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
*/


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

matrix_t* matrix_create_multiply(matrix_t* A, matrix_t* B) {
	assert(A != NULL);
	assert(B != NULL);
	assert(matrix_cols(A) == matrix_rows(B));

	matrix_t* output = matrix_create_constant(matrix_rows(A),matrix_cols(B),0.0);
	matrix_mutate_gemm(1, A, B, 0, output);

	return output;
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
/* END OF FILE                                                                */
/******************************************************************************/
