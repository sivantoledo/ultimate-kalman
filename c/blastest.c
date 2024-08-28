
/*
 * Sivan Toledo, 2024
 *
 * To generate the random numbers in Matlab, use
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifdef _WIN32
// for "unused" attribute
#define __attribute__(x)
#include <float.h>
// string.h for memcpy
#include <string.h>
#else
#include <unistd.h>
#endif

#include <math.h>

#include "ultimatekalman.h"

double A_rowwise[] = {
  1, 2, 3,
  4, 5, 6
};

double B_rowwise[] = {
  7,  8,
  9, 10,
 12, 13
};

int main(int argc, char* argv[]) {
	printf("BLAS test starting\n");

	kalman_matrix_t* C = matrix_create_constant(2, 2, 1.0);
	kalman_matrix_t* A = matrix_create_from_rowwise(A_rowwise, 2, 3);
	kalman_matrix_t* B = matrix_create_from_rowwise(B_rowwise, 3, 2);
	
	printf("A = ");
	matrix_print(A, "%.0f");

	printf("B = ");
	matrix_print(B, "%.0f");

	matrix_mutate_gemm(2.0, A, B, 3.0, C);

	printf("C = ");
	matrix_print(C, "%.0f");

	printf("Result should be:\n  125  137  \n  293  323\n");

	matrix_free(A);
	matrix_free(B);
	matrix_free(C);

	printf("BLAS test done\n");
	return 0;
}

     
