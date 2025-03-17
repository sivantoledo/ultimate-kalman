/*
 * kalman_explicit_representation.c
 *
 * (C) Sivan Toledo and Shahaf Gargir, 2024-2025
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// for getenv
#include <stdlib.h>

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
#include "kalman.h"
#include "parallel.h"
#include "concurrent_set.h"
#include "memory.h"

/******************************************************************************/
/* UTILITIES                                                                  */
/******************************************************************************/


/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

#ifdef moved
typedef struct step_st {
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
} step_t;
#endif

typedef kalman_step_equations_t step_t;

static void* step_create() {
  step_t *s = malloc(sizeof(step_t));
  s->step = -1;
  s->dimension = -1;

  s->C_type = 0;
  s->K_type = 0;

  s->H = NULL;
  s->F = NULL;
  s->K = NULL;
  s->c = NULL;

  s->G = NULL;
  s->o = NULL;
  s->C = NULL;

  s->state = NULL;
  s->covariance = NULL;
  s->covariance_type = 'C';

  assert(s != NULL);
  return s;
}

static void step_free(void *v) {
  step_t *s = (step_t*) v;

  matrix_free(s->H);
  matrix_free(s->F);
  matrix_free(s->K);
  matrix_free(s->c);

  matrix_free(s->G);
  matrix_free(s->o);
  matrix_free(s->C);

  matrix_free(s->state);
  matrix_free(s->covariance);

  free(s);
}

static void step_rollback(void *v) {
  step_t *s = (step_t*) v;

  matrix_free(s->G);
  matrix_free(s->o);
  matrix_free(s->C);

  matrix_free(s->state);
  matrix_free(s->covariance);
}

static kalman_step_index_t step_get_index(void *v) {
  return ((step_t*) v)->step;
}

static int32_t step_get_dimension(void *v) {
  return ((step_t*) v)->dimension;
}

static kalman_matrix_t* step_get_state(void *v) {
  return ((step_t*) v)->state;
}

static kalman_matrix_t* step_get_covariance(void *v) {
  return ((step_t*) v)->covariance;
}

static char step_get_covariance_type(void *v) {
  return ((step_t*) v)->covariance_type;
}

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

static void evolve(kalman_t *kalman, int32_t n_i, matrix_t *H_i, matrix_t *F_i, matrix_t *c_i, matrix_t *K_i, char K_type) {
  step_t *kalman_current;
  kalman->current = kalman_current = step_create();
  kalman_current->dimension = n_i;

  if (farray_size(kalman->steps) == 0) {
    //if (debug) printf("kalman_evolve first step\n");
    kalman_current->step = 0;
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

  step_t *imo = farray_get_last(kalman->steps);
  assert(imo != NULL);
  kalman_current->step = (imo->step) + 1;

  //if (debug) printf("kalman_evolve step = %d\n",kalman->current->step);

  assert(H_i != NULL);
  assert(F_i != NULL);
  assert(c_i != NULL);
  //assert(K_i!=NULL);

  // matrix_t* V_i_H_i = kalman_covariance_matrix_weigh(K_i,K_type,H_i);
  // matrix_t* V_i_F_i = kalman_covariance_matrix_weigh(K_i,K_type,F_i);
  // matrix_t* V_i_c_i = kalman_covariance_matrix_weigh(K_i,K_type,c_i);

  // matrix_mutate_scale(V_i_F_i,-1.0);

  kalman_current->H = matrix_create_copy(H_i);
  kalman_current->F = matrix_create_copy(F_i);
  kalman_current->c = matrix_create_copy(c_i);
  kalman_current->K = matrix_create_copy(K_i);
  kalman_current->K_type = K_type;
}

static void observe(kalman_t *kalman, matrix_t *G_i, matrix_t *o_i, matrix_t *C_i, char C_type) {

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe!\n");
#endif

  assert(kalman != NULL);
  assert(kalman->current != NULL);
  step_t *kalman_current = kalman->current;
  int32_t n_i = kalman_current->dimension;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe %d\n",(int) kalman->current->step);
#endif

  if (o_i != NULL) {
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("C_i(%c) ",C_type);
		matrix_print(C_i,"%.3e");
		printf("G_i ");
		matrix_print(G_i,"%.3e");
		printf("o_i ");
		matrix_print(o_i,"%.3e");
#endif
    kalman_current->G = matrix_create_copy(G_i);
    kalman_current->o = matrix_create_copy(o_i);
    kalman_current->C = matrix_create_copy(C_i);
    kalman_current->C_type = C_type;
  }

  farray_append(kalman->steps, kalman->current);
  kalman->current = NULL;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_observe done\n");
#endif
}

static void smooth(kalman_t *kalman) {
  fprintf(stderr,"explicit rep smooth\n");

  if (kalman->options & KALMAN_ALGORITHM_ODDEVEN)     kalman_smooth_oddeven    (kalman->options, (kalman_step_equations_t**) kalman->steps->elements, farray_size(kalman->steps));
  if (kalman->options & KALMAN_ALGORITHM_ASSOCIATIVE) kalman_smooth_associative(kalman->options, (kalman_step_equations_t**) kalman->steps->elements, farray_size(kalman->steps));
}

void kalman_create_explicit_representation(kalman_t *kalman) {
  kalman->evolve = evolve;
  kalman->observe = observe;
  kalman->smooth = smooth;

  kalman->step_create = step_create;
  kalman->step_free = step_free;
  kalman->step_rollback = step_rollback;
  kalman->step_get_index = step_get_index;
  kalman->step_get_dimension = step_get_dimension;
  kalman->step_get_state = step_get_state;
  kalman->step_get_covariance = step_get_covariance;
  kalman->step_get_covariance_type = step_get_covariance_type;
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
