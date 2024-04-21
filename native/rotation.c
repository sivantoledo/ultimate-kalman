
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

double PI = 3.141592653589793;

int main(int argc, char* argv[]) {

	printf("rotation starting\n");

	/*
	double alpha = 2.0 * PI / 16.0;

	double F_rowwise[] = {
			cos(alpha), -sin(alpha),
			sin(alpha),  cos(alpha),
	};

	kalman_matrix_t* F = matrix_create_from_rowwise(F_rowwise, 2, 2);

	matrix_print(F, "%.6f");
*/
	return 0;
}
