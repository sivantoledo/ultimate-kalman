/*
 * embarrasingly_parallel.c
 *
 * A program for testing the scalability of the building blocks of the
 * Gargir-Toledo parallel Kalman smoother, using parallel for loops
 * with no data dependencies. This is an ideal case whose speedups in
 * a given environment provide upper bounds on the speedups that the
 * Gargir-Toledo smoother is likely to achieve.
 *
 * Copyright (C) Sivan Toledo, 2024-2025
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// for getenv
//#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
int gettimeofday(struct timeval * tp, struct timezone * tzp);
#else
#include <sys/time.h>
#endif

#include "parallel.h"
#include "matrix_ops.h"
#include "cmdline_args.h"


/**************************************************************************************************/

double PI = 3.141592653589793;

// this function is from ChatGPT.
double generateGaussian(double mean, double stddev) {
  // Generate two uniform random numbers in the range (0, 1)
  double u1 = rand() / (RAND_MAX + 1.0);
  double u2 = rand() / (RAND_MAX + 1.0);

  // Apply Box-Muller transform
  double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

  // Adjust for desired mean and standard deviation
  return z0 * stddev + mean;
}

kalman_matrix_t* generateRandomOrthonormal(int rows, int cols) {
  int i, j;
  //int n;

  if (rows == 1 || cols == 1) {
    kalman_matrix_t *QR = matrix_create(rows, cols);
    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        matrix_set(QR, i, j, generateGaussian(0.0, 1.0));
      }
    }
    return QR;
  }

  // we assume that rows >= cols

  kalman_matrix_t *QR = matrix_create(rows, rows);
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      matrix_set(QR, i, j, generateGaussian(0.0, 1.0));
    }
  }

  kalman_matrix_t *A = matrix_create_identity(rows, cols);
  kalman_matrix_t *TAU = matrix_create_mutate_qr(QR);
  matrix_mutate_apply_qt(QR, TAU, A);

  matrix_free(TAU);
  matrix_free(QR);

  return A;
}

typedef struct step_st {
    kalman_matrix_t* A;
    kalman_matrix_t* TAU;
} step_t;

void ep_alloc_struct(void *indices_v, parallel_index_t n, parallel_index_t start, parallel_index_t end) {
  step_t** indices = (step_t**) indices_v;

  for (parallel_index_t i = start; i < end; ++i) {
    indices[i] = (step_t*) malloc(sizeof(step_t));
  }
}

void ep_alloc_matrix(void *indices_v, parallel_index_t n, parallel_index_t start, parallel_index_t end) {
  step_t **indices = (step_t**) indices_v;

  for (parallel_index_t i = start; i < end; ++i) {
    indices[i]->A = matrix_create((int32_t) (2 * n), (int32_t) n);
  }
}

void ep_fill_matrix(void *indices_v, parallel_index_t n, parallel_index_t start, parallel_index_t end) {
  step_t **indices = (step_t**) indices_v;

  for (parallel_index_t i = start; i < end; ++i) {
    for (int r = 0; r < 2 * n; r++) {
      for (int c = 0; c < n; c++) {
        matrix_set(indices[i]->A, r, c, (double) (r + c));
      }
    }
  }
}

void ep_create(void *indices_v, parallel_index_t n, parallel_index_t start, parallel_index_t end) {
  step_t **indices = (step_t**) indices_v;

  for (parallel_index_t i = start; i < end; ++i) {
    indices[i] = (step_t*) malloc(sizeof(step_t));
    //indices[i] -> A = matrix_create_copy(A);
    indices[i]->A = matrix_create((int32_t) (2 * n), (int32_t) n);
    for (int r = 0; r < 2 * n; r++) {
      for (int c = 0; c < n; c++) {
        matrix_set(indices[i]->A, r, c, (double) (r + c));
      }
    }
  }
}

void ep_factor(void *indices_v, parallel_index_t n, parallel_index_t start, parallel_index_t end) {
  step_t **indices = (step_t**) indices_v;

  for (parallel_index_t i = start; i < end; ++i) {
    indices[i]->TAU = matrix_create_mutate_qr(indices[i]->A);
  }
}

int main(int argc, char *argv[]) {

  double times[16];

  int n, k;
  int nthreads, blocksize;
  int present;

  parse_args(argc, argv);
  present = get_int_param    ("n",         &n, 6);
  present = get_int_param    ("k",         &k, 100000);
  present = get_int_param    ("nthreads",  &nthreads,  -1);
  present = get_int_param    ("blocksize", &blocksize, -1);
  check_unused_args();

  printf("embarrassingly_parallel n=dimension=%d k=count=%d nthreads=%d blocksize=%d (-1 means do not set)\n",n,k,nthreads,blocksize);

  if (nthreads != -1)  parallel_set_thread_limit(nthreads);
  if (blocksize != -1) parallel_set_blocksize(blocksize);

  double t = 0.0;

  struct timeval begin, end;
  long seconds, microseconds;

  gettimeofday(&begin, 0);

  step_t **indices = (step_t**) malloc(k * sizeof(step_t*));

  foreach_in_range(ep_alloc_struct, indices, n, k);

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  times[0] = seconds + microseconds * 1e-6;

  foreach_in_range(ep_alloc_matrix, indices, n, k);

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  times[1] = seconds + microseconds * 1e-6;

  foreach_in_range(ep_fill_matrix, indices, n, k);

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  times[2] = seconds + microseconds * 1e-6;

  foreach_in_range(ep_factor, indices, n, k);

  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  microseconds = end.tv_usec - begin.tv_usec;
  times[3] = seconds + microseconds * 1e-6;

  printf("embarrassingly parallel testing took %.2e seconds\n", times[3]);

  printf("performance testing breakdown %.2e %.2e %.2e %.2e (alloc_struct, alloc_matrix, fill_matrix, factor)\n",
      times[0], times[1] - times[0], times[2] - times[1], times[3] - times[2]);

  printf("embarrassingly parallel testing done\n");

  return 0;
}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
