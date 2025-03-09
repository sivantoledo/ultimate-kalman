#ifndef ULTIMATE_KALMAN_H
#define ULTIMATE_KALMAN_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef BUILD_MEX
#include "mex.h"
#endif

#include "kalman_matrix_ops.h"
#include "flexible_arrays.h"

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

//typedef struct cov_st {
	//int64_t first; // logical index of first element
//} cov_t;

//kalman_matrix_t* cov_weigh(kalman_matrix_t* cov, char cov_type, kalman_matrix_t* A);
//kalman_matrix_t* cov_weigh(cov_t* cov, kalman_matrix_t* A);

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

typedef struct step_st {
	int64_t step; // logical step number
	int32_t dimension;

	kalman_matrix_t* Rdiag;
	kalman_matrix_t* Rsupdiag;
	kalman_matrix_t* y;

	kalman_matrix_t* Rbar;
	kalman_matrix_t* ybar;

	kalman_matrix_t* state;
	kalman_matrix_t* covariance;
} step_t;

// Sivan Nov 2022: comment out since these are static functions.
//step_t* step_create();
//void    step_free(step_t* step);

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

void      kalman_evolve    (kalman_t* kalman, int32_t n_i, kalman_matrix_t* H_i, kalman_matrix_t* F_i, kalman_matrix_t* c_i, kalman_matrix_t* K_i, char K_type);
void      kalman_observe   (kalman_t* kalman, kalman_matrix_t* G_i, kalman_matrix_t* o_i, kalman_matrix_t* C_i, char C_type);
void      kalman_smooth    (kalman_t* kalman);
kalman_matrix_t* kalman_estimate  (kalman_t* kalman, int64_t si);
kalman_matrix_t* kalman_covariance(kalman_t* kalman, int64_t si);
char             kalman_covariance_type(kalman_t* kalman, int64_t si);
void      kalman_forget    (kalman_t* kalman, int64_t si);
void      kalman_rollback  (kalman_t* kalman, int64_t si);
kalman_matrix_t* kalman_perftest(kalman_t* kalman,
		                      kalman_matrix_t* H, kalman_matrix_t* F, kalman_matrix_t* c, kalman_matrix_t* K, char K_type,
		                      kalman_matrix_t* G, kalman_matrix_t* o,              kalman_matrix_t* C, char C_type,
									        int32_t count, int32_t decimation);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
