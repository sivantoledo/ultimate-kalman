/*
 * ultimatekalman.c
 *
 * (C) Sivan Toledo, 2022
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
#include "ultimatekalman.h"

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
		printf("unknown covariance-matrix type %c\n",cov_type);
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

typedef struct step_st {
	int64_t step; // logical step number
	int32_t dimension;
	
	char   C_type;
	char   K_type;

	//kalman_matrix_t* H;
	kalman_matrix_t* F;

	kalman_matrix_t* predictedState;
	kalman_matrix_t* predictedCovariance;

	kalman_matrix_t* assimilatedState;
	kalman_matrix_t* assimilatedCovariance;

	kalman_matrix_t* smoothedState;
	kalman_matrix_t* smoothedCovariance;

	kalman_matrix_t* state;
	kalman_matrix_t* covariance;
} step_t;

static step_t* step_create() {
	step_t* s = malloc(sizeof(step_t));
	s->step      = -1;
	s->dimension = -1;
	
	s->C_type    = 0;
	s->K_type    = 0;


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

	assert( s != NULL );
	return s;
}

void step_free(step_t* s) {
	
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

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

kalman_t* kalman_create() {
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
	step_t* kalman_current;
	kalman->current = kalman_current = step_create();
	kalman_current->dimension = n_i;

	if (farray_size(kalman->steps)==0) {
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

	step_t* imo = farray_get_last( kalman->steps );
	assert(imo != NULL);
	kalman_current->step = (imo->step) + 1;

	//if (debug) printf("kalman_evolve step = %d\n",kalman->current->step);

	assert(H_i!=NULL);
	assert(F_i!=NULL);
	assert(c_i!=NULL);

	// we assume H_i is an identity, need to check in the final code
	//kalman->current->H	= matrix_create_copy(H_i);
	kalman_current->F 	= matrix_create_copy(F_i);

	//printf("imo->assimilatedState = ");
	//matrix_print(imo->assimilatedState,"%.3e");
	matrix_t* predictedState = matrix_create_multiply( F_i, imo->assimilatedState );

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

	kalman_current->predictedState = matrix_create_add(predictedState,c_i);
	matrix_free(predictedState);

	matrix_t* K_i_explicit = explicit(K_i,K_type);

	matrix_t* t4 = matrix_create_multiply( F_i, imo->assimilatedCovariance );
	matrix_t* F_iTrans = matrix_create_transpose( F_i );
	matrix_t* t5 = matrix_create_multiply( t4, F_iTrans );

	kalman_current->predictedCovariance = matrix_create_add( t5, K_i_explicit );

	matrix_free( F_iTrans );
	matrix_free( t5 );
	matrix_free( t4 );

	kalman_current->state      = kalman_current->predictedState      ;
	kalman_current->covariance = kalman_current->predictedCovariance ;
}

void kalman_observe(kalman_t* kalman, matrix_t* G_i, matrix_t* o_i, matrix_t* C_i, char C_type) {

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe!\n");
#endif

	assert(kalman != NULL);
	assert(kalman->current != NULL);
	step_t* kalman_current = kalman->current;	
	int32_t n_i = kalman_current->dimension;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("observe %d\n",(int) kalman_current->step);
#endif

	// matrix_t* W_i_G_i = NULL;
	// matrix_t* W_i_o_i = NULL;

	if (kalman_current->step == 0) {
		//printf("filter_smoother step 0 observation cov-type %c\n",C_type);
		matrix_t* W_i_G_i = cov_weigh(C_i,C_type,G_i);
		matrix_t* W_i_o_i = cov_weigh(C_i,C_type,o_i);

		//matrix_print(W_i_G_i,"%.3e");
		//matrix_print(W_i_o_i,"%.3e");

		matrix_t* tau = matrix_create_mutate_qr( W_i_G_i );
		matrix_mutate_apply_qt( W_i_G_i, tau, W_i_o_i );

		matrix_mutate_chop(W_i_G_i, n_i, n_i); // it might have been tall
		matrix_mutate_chop(W_i_o_i, n_i, 1  ); // it might have been tall
		matrix_mutate_triu(W_i_G_i);

		kalman_current->assimilatedState =  matrix_create_trisolve(W_i_G_i,W_i_o_i);

		matrix_t* R_trans = matrix_create_transpose(W_i_G_i);
		matrix_t* R_trans_R = matrix_create_multiply( R_trans, W_i_G_i );

		kalman_current->assimilatedCovariance = matrix_create_inverse(R_trans_R);

		//matrix_print(kalman->current->assimilatedCovariance,"%.3e");
		//matrix_print(kalman->current->assimilatedState,"%.3e");

		matrix_free( R_trans_R );
		matrix_free( R_trans );
		matrix_free( tau );
		matrix_free( W_i_G_i );
		matrix_free( W_i_o_i );

		kalman_current->state      = kalman_current->assimilatedState      ;
		kalman_current->covariance = kalman_current->assimilatedCovariance ;

		farray_append(kalman->steps,kalman->current);
		return;
	}

	if (o_i == NULL) {
		kalman_current->assimilatedState      = matrix_create_copy( kalman_current->predictedState      );
		kalman_current->assimilatedCovariance = matrix_create_copy( kalman_current->predictedCovariance );

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

		matrix_t* predictedObservations = matrix_create_multiply( G_i, kalman_current->predictedState );

		matrix_t* G_i_trans = matrix_create_transpose( G_i );
		matrix_t* t1 = matrix_create_multiply( G_i, kalman_current->predictedCovariance );
		matrix_t* t2 = matrix_create_multiply( t1, G_i_trans );
		matrix_t* C_i_explicit = explicit( C_i, C_type );
		matrix_t* S = matrix_create_add( t2, C_i_explicit );

		matrix_t* t4 = matrix_create_multiply( kalman_current->predictedCovariance, G_i_trans );
		matrix_t* S_inv = matrix_create_inverse( S );
		matrix_t* gain = matrix_create_multiply( t4, S_inv );

		matrix_t* innovation = matrix_create_subtract( o_i, predictedObservations );

		matrix_t* gain_innovation = matrix_create_multiply( gain, innovation );

		kalman_current->assimilatedState = matrix_create_add( kalman_current->predictedState, gain_innovation );

		matrix_t* t5 = matrix_create_multiply( gain, G_i );
		matrix_t* t6 = matrix_create_multiply( t5, kalman_current->predictedCovariance );

		kalman_current->assimilatedCovariance = matrix_create_subtract( kalman_current->predictedCovariance, t6 );

		matrix_free( t6 );
		matrix_free( gain_innovation );
		matrix_free( innovation );
		matrix_free( gain );
		matrix_free( S_inv );
		matrix_free( t5 );
		matrix_free( t4 );
		matrix_free( S );
		matrix_free( C_i_explicit );
		//matrix_free( t3 );
		matrix_free( t2 );
		matrix_free( t1 );
		matrix_free( G_i_trans );
		matrix_free( predictedObservations );

	}

	kalman_current->state      = kalman_current->assimilatedState      ;
	kalman_current->covariance = kalman_current->assimilatedCovariance ;

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

/*
 * Sivan March 2025: this was called smooth and was called trivially from kalman_smooth
 */
void kalman_smooth(kalman_t* kalman) {
	if (farray_size(kalman->steps) == 0) return;

	int64_t si;
	int64_t last  = farray_last_index (kalman->steps);
	int64_t first = farray_first_index(kalman->steps);

	step_t* i = farray_get(kalman->steps,last);

	i->smoothedState      = matrix_create_copy( i->assimilatedState      );
	i->smoothedCovariance = matrix_create_copy( i->assimilatedCovariance );

	i->state      = i->smoothedState      ;
	i->covariance = i->smoothedCovariance ;

	printf("smooth first:last = %lld:%lld\n",first,last);
	step_t* ipo = i;
	for (si=last-1; si>=first; si--) {
		i = farray_get(kalman->steps,si);

		matrix_t* nextPredictedEstimate    = ipo->predictedState;
		matrix_t* nextPredictedCovariance  = ipo->predictedCovariance;

		matrix_t* nextPredictedCovarianceInv  = matrix_create_inverse( nextPredictedCovariance );

		matrix_t* nextSmoothEstimate       = ipo->smoothedState;
		matrix_t* nextSmoothCovariance     = ipo->smoothedCovariance;

		//matrix_t* nextEvolutionMatrix      = matrix_create_mldivide( ipo->H, ipo->F );
		matrix_t* nextEvolutionMatrix      = ipo->F;
		matrix_t* nextEvolutionMatrixTrans = matrix_create_transpose( nextEvolutionMatrix );

		matrix_t* assimilatedState         = i->assimilatedState;
		matrix_t* assimilatedCovariance    = i->assimilatedCovariance;

		matrix_t* t1                       = matrix_create_multiply( assimilatedCovariance, nextEvolutionMatrixTrans );
		matrix_t* backwardInnovation       = matrix_create_multiply( t1, nextPredictedCovarianceInv );
		matrix_t* backwardInnovationTrans  = matrix_create_transpose( backwardInnovation );
		matrix_t* t2                       = matrix_create_subtract( nextSmoothEstimate, nextPredictedEstimate );
		matrix_t* t3                       = matrix_create_multiply( backwardInnovation, t2 );
		i->smoothedState                   = matrix_create_add( assimilatedState, t3 );

		matrix_t* t4                       = matrix_create_subtract( nextSmoothCovariance, nextPredictedCovariance );
		matrix_t* t5                       = matrix_create_multiply( backwardInnovation, t4 );
		matrix_t* t6                       = matrix_create_multiply( t5, backwardInnovationTrans );
		i->smoothedCovariance              = matrix_create_add( assimilatedCovariance, t6 );

		matrix_free( t6 );
		matrix_free( t5 );
		matrix_free( t4 );
		matrix_free( t3 );
		matrix_free( t2 );
		matrix_free( t1 );
		//matrix_free( nextEvolutionMatrix );
		matrix_free( nextEvolutionMatrixTrans );
		matrix_free( backwardInnovationTrans );
		matrix_free( backwardInnovation );
		matrix_free( nextPredictedCovarianceInv );

		i->state      = i->smoothedState      ;
		i->covariance = i->smoothedCovariance ;

		ipo = i;
	}

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

			matrix_free(step->smoothedState);
			matrix_free(step->smoothedCovariance);
			matrix_free(step->assimilatedState);
			matrix_free(step->assimilatedCovariance);
			// fix aliases
			step->state      = step->predictedState;
			step->covariance = step->predictedCovariance;

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
/* END OF FILE                                                                */
/******************************************************************************/
