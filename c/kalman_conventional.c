/*
 * kalman_filter_smoother.c
 *
 * (C) Sivan Toledo, 2024-2025
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

#define KALMAN_MATRIX_SHORT_TYPE
#include "kalman.h"

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

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

typedef struct step_st {
  int64_t step; // logical step number
  int32_t dimension;

  char C_type;
  char K_type;

  //kalman_matrix_t* H;
  kalman_matrix_t *F;

  kalman_matrix_t *predictedState;
  kalman_matrix_t *predictedCovariance;

  kalman_matrix_t *assimilatedState;
  kalman_matrix_t *assimilatedCovariance;

  kalman_matrix_t *smoothedState;
  kalman_matrix_t *smoothedCovariance;

  kalman_matrix_t *state;
  kalman_matrix_t *covariance;
} step_t;

//static step_t* step_create() {
static void* step_create() {
  step_t *s = malloc(sizeof(step_t));
  s->step = -1;
  s->dimension = -1;

  s->C_type = 0;
  s->K_type = 0;

  //s->H = NULL;
  s->F = NULL;

  s->predictedState = NULL;
  s->predictedCovariance = NULL;
  s->assimilatedState = NULL;
  s->assimilatedCovariance = NULL;
  s->smoothedState = NULL;
  s->smoothedCovariance = NULL;
  s->state = NULL;
  s->covariance = NULL;

  assert(s != NULL);
  return s;
}

//void step_free(step_t* s) {
static void step_free(void *v) {
  step_t *s = (step_t*) v;

  //matrix_free(s->H);
  matrix_free(s->F);

  matrix_free(s->predictedState);
  matrix_free(s->predictedCovariance);
  matrix_free(s->assimilatedState);
  matrix_free(s->assimilatedCovariance);
  matrix_free(s->smoothedState);
  matrix_free(s->smoothedCovariance);

  // state and covariance are aliases in this implementation

  free(s);
}

static void step_rollback(void *v) {
  step_t *s = (step_t*) v;
  matrix_free(s->smoothedState);
  matrix_free(s->smoothedCovariance);
  matrix_free(s->assimilatedState);
  matrix_free(s->assimilatedCovariance);
  // fix aliases
  s->state = s->predictedState;
  s->covariance = s->predictedCovariance;
}

static int64_t step_get_index(void *v) {
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
  return 'C';
}

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

static void evolve(kalman_t *kalman, int32_t n_i, matrix_t *H_i, matrix_t *F_i, matrix_t *c_i, matrix_t *K_i,
    char K_type) {
  step_t *kalman_current;
  kalman->current = kalman_current = (step_t*) step_create();
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

  // we assume H_i is an identity, need to check in the final code
  //kalman->current->H	= matrix_create_copy(H_i);
  kalman_current->F = matrix_create_copy(F_i);

  //printf("imo->assimilatedState = ");
  //matrix_print(imo->assimilatedState,"%.3e");
  matrix_t *predictedState = matrix_create_multiply(F_i, imo->assimilatedState);

  // with H_i
  //matrix_t* Hinv = matrix_create_inverse(H_i);
  //matrix_t* HinvTrans = matrix_create_transpose( Hinv );
  //matrix_t* predictedPlusCi = matrix_create_add(predictedState,c_i);
  //matrix_free(predictedState);
  //kalman->current->predictedState = matrix_create_multiply(Hinv, predictedPlusCi);
  //matrix_free(predictedPlusCi);

  //matrix_t* K_i_explicit = explicit(K_i,K_type);
  //matrix_t* t1 = matrix_create_multiply(Hinv,K_i_explicit);
  //matrix_t* t2 = matrix_create_multiply(t1,HinvTrans);

  //matrix_t* t3 = matrix_create_multiply( Hinv, F_i );
  //matrix_t* t4 = matrix_create_multiply( t3, imo->assimilatedCovariance );
  //matrix_t* F_iTrans = matrix_create_transpose( F_i );
  //matrix_t* t5 = matrix_create_multiply( t4, F_iTrans );
  //matrix_t* t6 = matrix_create_multiply( t5, HinvTrans );

  //matrix_free( t6 );
  //matrix_free( F_iTrans );
  //matrix_free( t5 );
  //matrix_free( t4 );
  //matrix_free( t3 );
  //matrix_free( t2 );
  //matrix_free( t1 );
  //matrix_free( HinvTrans );
  //matrix_free( Hinv );

  kalman_current->predictedState = matrix_create_add(predictedState, c_i);
  matrix_free(predictedState);

  matrix_t *K_i_explicit = explicit(K_i, K_type);

  matrix_t *t4 = matrix_create_multiply(F_i, imo->assimilatedCovariance);
  matrix_t *F_iTrans = matrix_create_transpose(F_i);
  matrix_t *t5 = matrix_create_multiply(t4, F_iTrans);

  kalman_current->predictedCovariance = matrix_create_add(t5, K_i_explicit);

  matrix_free(F_iTrans);
  matrix_free(t5);
  matrix_free(t4);

  kalman_current->state = kalman_current->predictedState;
  kalman_current->covariance = kalman_current->predictedCovariance;
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
	printf("observe %d\n",(int) kalman_current->step);
#endif

  // matrix_t* W_i_G_i = NULL;
  // matrix_t* W_i_o_i = NULL;

  if (kalman_current->step == 0) {
    //printf("filter_smoother step 0 observation cov-type %c\n",C_type);
    matrix_t *W_i_G_i = cov_weigh(C_i, C_type, G_i);
    matrix_t *W_i_o_i = cov_weigh(C_i, C_type, o_i);

    //matrix_print(W_i_G_i,"%.3e");
    //matrix_print(W_i_o_i,"%.3e");

    matrix_t *tau = matrix_create_mutate_qr(W_i_G_i);
    matrix_mutate_apply_qt(W_i_G_i, tau, W_i_o_i);

    matrix_mutate_chop(W_i_G_i, n_i, n_i); // it might have been tall
    matrix_mutate_chop(W_i_o_i, n_i, 1); // it might have been tall
    matrix_mutate_triu(W_i_G_i);

    kalman_current->assimilatedState = matrix_create_trisolve(W_i_G_i, W_i_o_i);

    matrix_t *R_trans = matrix_create_transpose(W_i_G_i);
    matrix_t *R_trans_R = matrix_create_multiply(R_trans, W_i_G_i);

    kalman_current->assimilatedCovariance = matrix_create_inverse(R_trans_R);

    //matrix_print(kalman->current->assimilatedCovariance,"%.3e");
    //matrix_print(kalman->current->assimilatedState,"%.3e");

    matrix_free(R_trans_R);
    matrix_free(R_trans);
    matrix_free(tau);
    matrix_free(W_i_G_i);
    matrix_free(W_i_o_i);

    kalman_current->state = kalman_current->assimilatedState;
    kalman_current->covariance = kalman_current->assimilatedCovariance;

    farray_append(kalman->steps, kalman->current);
    return;
  }

  if (o_i == NULL) {
    kalman_current->assimilatedState = matrix_create_copy(kalman_current->predictedState);
    kalman_current->assimilatedCovariance = matrix_create_copy(kalman_current->predictedCovariance);

  } else {
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

    matrix_t *predictedObservations = matrix_create_multiply(G_i, kalman_current->predictedState);

    matrix_t *G_i_trans = matrix_create_transpose(G_i);
    matrix_t *t1 = matrix_create_multiply(G_i, kalman_current->predictedCovariance);
    matrix_t *t2 = matrix_create_multiply(t1, G_i_trans);
    matrix_t *C_i_explicit = explicit(C_i, C_type);
    matrix_t *S = matrix_create_add(t2, C_i_explicit);

    matrix_t *t4 = matrix_create_multiply(kalman_current->predictedCovariance, G_i_trans);
    matrix_t *S_inv = matrix_create_inverse(S);
    matrix_t *gain = matrix_create_multiply(t4, S_inv);

    matrix_t *innovation = matrix_create_subtract(o_i, predictedObservations);

    matrix_t *gain_innovation = matrix_create_multiply(gain, innovation);

    kalman_current->assimilatedState = matrix_create_add(kalman_current->predictedState, gain_innovation);

    matrix_t *t5 = matrix_create_multiply(gain, G_i);
    matrix_t *t6 = matrix_create_multiply(t5, kalman_current->predictedCovariance);

    kalman_current->assimilatedCovariance = matrix_create_subtract(kalman_current->predictedCovariance, t6);

    matrix_free(t6);
    matrix_free(gain_innovation);
    matrix_free(innovation);
    matrix_free(gain);
    matrix_free(S_inv);
    matrix_free(t5);
    matrix_free(t4);
    matrix_free(S);
    matrix_free(C_i_explicit);
    //matrix_free( t3 );
    matrix_free(t2);
    matrix_free(t1);
    matrix_free(G_i_trans);
    matrix_free(predictedObservations);

  }

  kalman_current->state = kalman_current->assimilatedState;
  kalman_current->covariance = kalman_current->assimilatedCovariance;

  farray_append(kalman->steps, kalman->current);

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_observe done\n");
#endif
}

static void smooth(kalman_t *kalman) {
  if (farray_size(kalman->steps) == 0)
    return;

  int64_t si;
  int64_t last = farray_last_index(kalman->steps);
  int64_t first = farray_first_index(kalman->steps);

  step_t *i = farray_get(kalman->steps, last);

  i->smoothedState = matrix_create_copy(i->assimilatedState);
  i->smoothedCovariance = matrix_create_copy(i->assimilatedCovariance);

  i->state = i->smoothedState;
  i->covariance = i->smoothedCovariance;

  printf("smooth first:last = %lld:%lld\n", first, last);
  step_t *ipo = i;
  for (si = last - 1; si >= first; si--) {
    i = farray_get(kalman->steps, si);

    matrix_t *nextPredictedEstimate = ipo->predictedState;
    matrix_t *nextPredictedCovariance = ipo->predictedCovariance;

    matrix_t *nextPredictedCovarianceInv = matrix_create_inverse(nextPredictedCovariance);

    matrix_t *nextSmoothEstimate = ipo->smoothedState;
    matrix_t *nextSmoothCovariance = ipo->smoothedCovariance;

    //matrix_t* nextEvolutionMatrix      = matrix_create_mldivide( ipo->H, ipo->F );
    matrix_t *nextEvolutionMatrix = ipo->F;
    matrix_t *nextEvolutionMatrixTrans = matrix_create_transpose(nextEvolutionMatrix);

    matrix_t *assimilatedState = i->assimilatedState;
    matrix_t *assimilatedCovariance = i->assimilatedCovariance;

    matrix_t *t1 = matrix_create_multiply(assimilatedCovariance, nextEvolutionMatrixTrans);
    matrix_t *backwardInnovation = matrix_create_multiply(t1, nextPredictedCovarianceInv);
    matrix_t *backwardInnovationTrans = matrix_create_transpose(backwardInnovation);
    matrix_t *t2 = matrix_create_subtract(nextSmoothEstimate, nextPredictedEstimate);
    matrix_t *t3 = matrix_create_multiply(backwardInnovation, t2);
    i->smoothedState = matrix_create_add(assimilatedState, t3);

    matrix_t *t4 = matrix_create_subtract(nextSmoothCovariance, nextPredictedCovariance);
    matrix_t *t5 = matrix_create_multiply(backwardInnovation, t4);
    matrix_t *t6 = matrix_create_multiply(t5, backwardInnovationTrans);
    i->smoothedCovariance = matrix_create_add(assimilatedCovariance, t6);

    matrix_free(t6);
    matrix_free(t5);
    matrix_free(t4);
    matrix_free(t3);
    matrix_free(t2);
    matrix_free(t1);
    //matrix_free( nextEvolutionMatrix );
    matrix_free(nextEvolutionMatrixTrans);
    matrix_free(backwardInnovationTrans);
    matrix_free(backwardInnovation);
    matrix_free(nextPredictedCovarianceInv);

    i->state = i->smoothedState;
    i->covariance = i->smoothedCovariance;

    ipo = i;
  }

}

void kalman_create_conventional(kalman_t *kalman) {
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
