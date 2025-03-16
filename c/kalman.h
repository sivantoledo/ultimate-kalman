#ifndef ULTIMATE_KALMAN_H
#define ULTIMATE_KALMAN_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef BUILD_MEX
#include "mex.h"
#endif

extern double kalman_nan;
#include "matrix_ops.h"
#include "flexible_arrays.h"

#ifdef KALMAN_STEP_INDEX_TYPE_INT32
typedef int32_t kalman_step_index_t;
#endif

#ifdef KALMAN_STEP_INDEX_TYPE_UINT32
typedef uint32_t kalman_step_index_t;
#endif

#ifdef KALMAN_STEP_INDEX_TYPE_INT64
typedef int64_t kalman_step_index_t;
#endif

#ifdef KALMAN_STEP_INDEX_TYPE_UINT64
typedef uint64_t kalman_step_index_t;
#endif

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

kalman_matrix_t* kalman_covariance_matrix_weigh   (kalman_matrix_t *cov, char cov_type, kalman_matrix_t *A);
kalman_matrix_t* kalman_covariance_matrix_explicit(kalman_matrix_t* cov, char type);

/******************************************************************************/
/* STEPS                                                                      */
/******************************************************************************/

/*
 * This structure is used as the input and output for parallel smoothers.
 *
 * It is meant to be used directly by clients of parallel smoothers.
 */
typedef struct kalman_step_equations_st {
  kalman_step_index_t step; // logical step number
  int32_t dimension;

  char C_type;
  char K_type;

  kalman_matrix_t *H;
  kalman_matrix_t *F;
  kalman_matrix_t *K;
  kalman_matrix_t *c;

  kalman_matrix_t *G;
  kalman_matrix_t *o;
  kalman_matrix_t *C;

  kalman_matrix_t *state;
  kalman_matrix_t *covariance;
} kalman_step_equations_t;

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

typedef enum {
  KALMAN_NONE                      = 0,      // No flags
  KALMAN_ALGORITHM_ULTIMATE        = 1 << 0, // 0x01 (1)
  KALMAN_ALGORITHM_CONVENTIONAL    = 1 << 1, // 0x02 (2)
  KALMAN_ALGORITHM_ODDEVEN         = 1 << 2,
  KALMAN_ALGORITHM_ASSOCIATIVE     = 1 << 3,
  KALMAN_NO_COVARIANCE             = 1 << 16
} kalman_options_t;

struct kalman_st;

typedef struct kalman_st {
    farray_t *steps;
    void *current; // really a pointer to step_t, but step_t varies among implementations
    kalman_options_t options;

    // implementation-specific operations
    void (*evolve)(struct kalman_st *kalman, int32_t n_i, kalman_matrix_t *H_i, kalman_matrix_t *F_i,
        kalman_matrix_t *c_i, kalman_matrix_t *K_i, char K_type);
    void (*observe)(struct kalman_st *kalman, kalman_matrix_t *G_i, kalman_matrix_t *o_i, kalman_matrix_t *C_i,
        char C_type);
    void (*smooth)(struct kalman_st *kalman);

    // functions on steps
    void* (*step_create)();         // returns a pointer to a step_t
    void (*step_free)(void *step_v);  // takes a pointer to a step_t
    void (*step_rollback)(void *step_v);  // rollback the step to just after the call to evolve

    kalman_step_index_t (*step_get_index)(void *step_v);
    int32_t (*step_get_dimension)(void *step_v);
    kalman_matrix_t* (*step_get_state)(void *step_v);
    kalman_matrix_t* (*step_get_covariance)(void *step_v);
    char (*step_get_covariance_type)();

} kalman_t;

kalman_t* kalman_create();
kalman_t* kalman_create_options(uint32_t options);
void kalman_free(kalman_t *kalman);

kalman_step_index_t kalman_earliest(kalman_t *kalman);
kalman_step_index_t kalman_latest(kalman_t *kalman);

void kalman_evolve(kalman_t *kalman, int32_t n_i, kalman_matrix_t *H_i, kalman_matrix_t *F_i, kalman_matrix_t *c_i,
    kalman_matrix_t *K_i, char K_type);
void kalman_observe(kalman_t *kalman, kalman_matrix_t *G_i, kalman_matrix_t *o_i, kalman_matrix_t *C_i, char C_type);
void kalman_smooth(kalman_t *kalman);

kalman_matrix_t* kalman_estimate(kalman_t *kalman, kalman_step_index_t si);
kalman_matrix_t* kalman_covariance(kalman_t *kalman, kalman_step_index_t si);
char kalman_covariance_type(kalman_t *kalman, kalman_step_index_t si);
void kalman_forget(kalman_t *kalman, kalman_step_index_t si);
void kalman_rollback(kalman_t *kalman, kalman_step_index_t si);
kalman_matrix_t* kalman_perftest(kalman_t *kalman,
                                 kalman_matrix_t *H,
                                 kalman_matrix_t *F,
                                 kalman_matrix_t *c,
                                 kalman_matrix_t *K, char K_type,
                                 kalman_matrix_t *G,
                                 kalman_matrix_t *o,
                                 kalman_matrix_t *C, char C_type,
                                 int32_t count, int32_t decimation);

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/

#endif
