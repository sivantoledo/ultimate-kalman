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

#ifdef _WIN32
#include <windows.h>
int gettimeofday(struct timeval * tp, struct timezone * tzp);
#else
#include <sys/time.h>
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
#include "memory.h"

double kalman_nan = 0.0 / 0.0;

/******************************************************************************/
/* COVARIANCE MATRICES                                                        */
/******************************************************************************/

matrix_t* kalman_covariance_matrix_weigh(matrix_t *cov, char cov_type, matrix_t *A) {
  blas_int_t M, N, K, LDA, LDB, LDC, NRHS, INFO;
  double ALPHA, BETA;

  assert(A != NULL);
  assert(cov != NULL);

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("cov(%c) ",cov_type);
	matrix_print(cov,"%.3e");
	printf("A ");
	matrix_print(A,"%.3e");
#endif

  matrix_t *WA = NULL;

  int32_t i, j, rows, cols;

  switch (cov_type) {
    case 'W':
      //if (debug) printf("cov W %d %d %d %d\n",matrix_cols(cov),matrix_rows(cov),matrix_cols(A),matrix_rows(A));

      assert(matrix_cols(cov) == matrix_rows(A));

      WA = matrix_create_constant(matrix_rows(A), matrix_cols(A), 0.0);

      M = matrix_rows(WA);
      N = matrix_cols(WA);
      K = matrix_cols(cov);
      LDA = matrix_ld(cov);
      LDB = matrix_ld(A);
      LDC = matrix_ld(WA);
      ALPHA = 1.0;
      BETA = 0.0;
#ifdef BUILD_BLAS_UNDERSCORE
    dgemm_
#else
      dgemm
#endif
      ("n", "N", &M, &N, &K, &ALPHA, cov->elements, &LDA, A->elements, &LDB, &BETA, WA->elements, &LDC
#ifdef BUILD_BLAS_STRLEN_END
          ,1,1
#endif
          );
      break;
    case 'U': // cov and an upper triangular matrix that we need to solve with
    case 'F': // same representation; in Matlab, we started from explicit cov and factored

      //if (debug) printf("cov U %d %d %d %d\n",matrix_cols(cov),matrix_rows(cov),matrix_cols(A),matrix_rows(A));

      // is this correct for 'F'?
      WA = matrix_create_trisolve(cov, A);
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

      WA = matrix_create_constant(rows, cols, 0.0);

      for (j = 0; j < cols; j++) {
        for (i = 0; i < rows; i++) {
          matrix_set(WA, i, j, matrix_get(cov, i, 0) * matrix_get(A, i, j));
        }
      }
      break;
    default:
      printf("unknown covariance-matrix type %c\n", cov_type);
      assert(0);
      WA = matrix_create_constant(matrix_rows(A), matrix_cols(A), kalman_nan);
      break;
  }

  //if (debug) printf("WA ");
  //if (debug) matrix_print(WA,"%.3e");

  return WA;
}

// SUPPORT 'W','C' only
matrix_t* kalman_covariance_matrix_explicit(matrix_t *cov, char type) {
  //fprintf(stderr,"cov type %c\n",type);
  assert(type == 'w' || type == 'W' || type == 'C' || type == 'U' || type == 'F');
  if (type == 'w') {
    blas_int_t n = matrix_rows(cov);

    matrix_t* C = matrix_create_constant(n, n, 0.0);

    for (int j = 0; j < n; j++) {
      matrix_set(C, j, j, pow(matrix_get(cov, j, 0),-2.0) );
    }
    return C;
  }
  if (type == 'W') {
    matrix_t *WT = matrix_create_transpose(cov);
    matrix_t *WTW = matrix_create_multiply(WT, cov);
    matrix_t *C = matrix_create_inverse(WTW);
    matrix_free(WT);
    matrix_free(WTW);
    return C;
  }
  if (type == 'U' || type == 'F') {
    //printf("L = ");
    //matrix_print(cov,"%.3e");
    matrix_t *LT = matrix_create_transpose(cov);
    matrix_t *C = matrix_create_multiply(cov, LT);
    matrix_free(LT);
    //printf("C = ");
    //matrix_print(C,"%.3e");
    return C;
  }
  if (type == 'C') {
    matrix_t *C = matrix_create_copy(cov);
    return C;
  }
  return NULL;
}

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

/* algorithm-specific constructors */
void kalman_create_ultimate    (kalman_t*);
void kalman_create_conventional(kalman_t*);
void kalman_create_oddeven     (kalman_t*);
void kalman_create_associative (kalman_t*);
void kalman_create_explicit_representation(kalman_t*);

kalman_t* kalman_create() {
  return kalman_create_options( KALMAN_ALGORITHM_ULTIMATE ); // default
  //return kalman_create_options( KALMAN_ALGORITHM_ULTIMATE | KALMAN_NO_COVARIANCE ); // default
  //return kalman_create_options( KALMAN_ALGORITHM_ODDEVEN | KALMAN_NO_COVARIANCE ); // default
  //return kalman_create_options(KALMAN_ALGORITHM_ASSOCIATIVE); // default
}

kalman_t* kalman_create_options(kalman_options_t options) {
#ifdef MOVED_TO_CLIENTS
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

  kalman_t *kalman = malloc(sizeof(kalman_t));
  assert(kalman != NULL);
  kalman->steps = farray_create();
  kalman->current = NULL;
  kalman->options = options;

  switch (options & (KALMAN_ALGORITHM_ULTIMATE | KALMAN_ALGORITHM_CONVENTIONAL | KALMAN_ALGORITHM_ODDEVEN | KALMAN_ALGORITHM_ASSOCIATIVE)) {
    case KALMAN_ALGORITHM_ULTIMATE:
      printf("calling kalman_ultimate\n");
      kalman_create_ultimate(kalman);
      break;
    case KALMAN_ALGORITHM_CONVENTIONAL:
      printf("calling kalman_filter_smoother\n");
      kalman_create_conventional(kalman);
      break;
    case KALMAN_ALGORITHM_ODDEVEN:
      printf("calling kalman_oddeven\n");
      //kalman_create_oddeven(kalman);
      kalman_create_explicit_representation(kalman);
      break;
    case KALMAN_ALGORITHM_ASSOCIATIVE:
      printf("calling kalman_associative\n");
      //kalman_create_associative(kalman);
      kalman_create_explicit_representation(kalman);
      break;
  }

  return kalman;
}

void kalman_free(kalman_t *kalman) {
  //printf("waning: kalman_free not fully implemented yet (steps not processed)\n");

  while (farray_size(kalman->steps) > 0) {
    void *i = farray_drop_last(kalman->steps);
    (*(kalman->step_free))(i);
  }

  farray_free(kalman->steps);
  // step_free( kalman->current );
  free(kalman);
}

kalman_step_index_t kalman_earliest(kalman_t *kalman) {
  if (farray_size(kalman->steps) == 0)
    return -1;
  //step_t* s = farray_get_first(kalman->steps);
  //return s->step;
  void *s = farray_get_first(kalman->steps);
  return (*(kalman->step_get_index))(s);
}

kalman_step_index_t kalman_latest(kalman_t *kalman) {
  if (farray_size(kalman->steps) == 0)
    return -1;
  //step_t* s = farray_get_last(kalman->steps);
  //return s->step;
  void *s = farray_get_last(kalman->steps);
  return (*(kalman->step_get_index))(s);
}

void kalman_evolve(kalman_t *kalman, int32_t n_i, matrix_t *H_i, matrix_t *F_i, matrix_t *c_i, matrix_t *K_i,
    char K_type) {
  (*(kalman->evolve))(kalman, n_i, H_i, F_i, c_i, K_i, K_type);
}

void kalman_observe(kalman_t *kalman, matrix_t *G_i, matrix_t *o_i, matrix_t *C_i, char C_type) {
  (*(kalman->observe))(kalman, G_i, o_i, C_i, C_type);
}

void kalman_smooth(kalman_t *kalman) {
  (*(kalman->smooth))(kalman);
}

matrix_t* kalman_estimate(kalman_t *kalman, kalman_step_index_t si) {
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_estimate\n");
#endif

  if (farray_size(kalman->steps) == 0)
    return NULL;

  if (si < 0)
    si = farray_last_index(kalman->steps);
  void *step = farray_get(kalman->steps, si);
  matrix_t *state = (*(kalman->step_get_state))(step);

  if (state == NULL) {
    int32_t dim = (*(kalman->step_get_dimension))(step);
    return matrix_create_constant(dim, 1, kalman_nan);
  }

  return matrix_create_copy(state);
}

void kalman_forget(kalman_t *kalman, kalman_step_index_t si) {
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d\n",(int) si);
#endif

  if (farray_size(kalman->steps) == 0)
    return;

  if (si < 0)
    si = farray_last_index(kalman->steps) - 1;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d now %d to %d\n",(int) si,(int) farray_first_index(kalman->steps),(int) farray_last_index(kalman->steps));
#endif

  if (si > farray_last_index(kalman->steps) - 1)
    return; // don't delete last step
  if (si < farray_first_index(kalman->steps))
    return; // nothing to delete

  while (farray_first_index(kalman->steps) <= si) {
    void *step = farray_drop_first(kalman->steps);
    (*(kalman->step_free))(step);
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("forget new first %d\n",(int) farray_first_index(kalman->steps));
#endif
  }
#ifdef BUILD_DEBUG_PRINTOUTS
	printf("forget %d done\n",(int) si);
#endif
}

void kalman_rollback(kalman_t *kalman, kalman_step_index_t si) {
  //printf("rollback %d\n",si);
  if (farray_size(kalman->steps) == 0)
    return;

  if (si > farray_last_index(kalman->steps))
    return; // we can roll  back even the last step (its observation)
  if (si < farray_first_index(kalman->steps))
    return;

  //step_t* step;
  void *step;
  kalman_step_index_t sj;
  do {
    step = farray_drop_last(kalman->steps);
    sj = (*(kalman->step_get_index))(step);
    if (sj == si) {
      //printf("rollback got to %d\n",si);

      /*
       matrix_free( step->Rdiag     );
       matrix_free( step->Rsupdiag  );
       matrix_free( step->y         );
       matrix_free( step->state     );
       matrix_free( step->covariance);
       */

      (*(kalman->step_rollback))(step);

      kalman->current = step;

    } else {
      //printf("rollback dropped step %d\n",step_get_step(step));
      // don't we need to free the step? Sivan March 2025
      (*(kalman->step_free))(step);
    }
  } while (sj > si);

  //printf("rollback to %d new latest %d\n",si,kalman_latest(kalman));

}

char kalman_covariance_type(kalman_t *kalman, kalman_step_index_t si) {
  //return (*(kalman->step_get_covariance_type))(); // currently the same for all steps

  if (farray_size(kalman->steps) == 0)
    return NULL;

  if (si < 0)
    si = farray_last_index(kalman->steps);
  void *step = farray_get(kalman->steps, si);

  return (*(kalman->step_get_covariance_type))(step);
}

matrix_t* kalman_covariance(kalman_t *kalman, kalman_step_index_t si) {

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_covariance\n");
#endif

  if (farray_size(kalman->steps) == 0)
    return NULL;

  if (si < 0)
    si = farray_last_index(kalman->steps);
  void *step = farray_get(kalman->steps, si);

  matrix_t *cov = NULL;

  if ((*(kalman->step_get_covariance))(step) != NULL) {
    cov = matrix_create_copy((*(kalman->step_get_covariance))(step));
  } else {
    int32_t n_i = (*(kalman->step_get_dimension))(step);
    cov = matrix_create_constant(n_i, n_i, kalman_nan);
  }

  return cov;
}

static struct timeval begin, end;

matrix_t* kalman_perftest(kalman_t *kalman, matrix_t *H, matrix_t *F, matrix_t *c, matrix_t *K, char K_type,
    matrix_t *G, matrix_t *o, matrix_t *C, char C_type, int32_t count, int32_t decimation) {

  matrix_t *t = matrix_create_constant(count / decimation, 1, kalman_nan);
  int32_t i, j, n;

  //printf("perftest count %d decimation %d rows %d\n",count,decimation,matrix_rows(t));

  //struct timeval begin, end;
  gettimeofday(&begin, 0);

  j = 0;
  n = matrix_cols(G);

  for (i = 0; i < count; i++) {
    //printf("perftest iter %d (j=%d)\n",i,j);
    //if (debug) printf("perftest iter %d (j=%d)\n",i,j);
    kalman_evolve(kalman, n, H, F, c, K, K_type);
    kalman_observe(kalman, G, o, C, C_type);
    matrix_t *e = kalman_estimate(kalman, -1);
    matrix_free(e);
    kalman_forget(kalman, -1);

    if ((i % decimation) == decimation - 1) {
      gettimeofday(&end, 0);
      long seconds = end.tv_sec - begin.tv_sec;
      long microseconds = end.tv_usec - begin.tv_usec;
      double elapsed = seconds + microseconds * 1e-6;

      assert(j < matrix_rows(t));

      //printf("perftest store i=%d (j=%d) %.3e\n",i,j,elapsed/decimation);
      matrix_set(t, j, 0, elapsed / decimation); // average per step
      j = j + 1;

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
  	printf("farray get(%d) %lld\n",i,(kalman_step_index_t) farray_get(a,i));
  }

  farray_drop_first(a);
  farray_drop_first(a);

  printf("farray size %lld\n",farray_size(a));
  printf("farray first %lld\n",farray_first_index(a));
  printf("farray last  %lld\n",farray_last_index(a));
  for (int i=farray_first_index(a); i<=farray_last_index(a); i++) {
  	printf("farray get(%d) %lld\n",i,(kalman_step_index_t) farray_get(a,i));
  }

	printf("farray get_first %lld\n",(kalman_step_index_t) farray_get_first(a));
	printf("farray get_last  %lld\n",(kalman_step_index_t) farray_get_last(a));

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
