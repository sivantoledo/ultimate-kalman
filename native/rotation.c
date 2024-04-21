
/*
 * Sivan Toledo, 2024
 *
 * To generate the random numbers in Matlab, use
 * rng(5); for j=2:16; evolErrs(1:2,j-1) = randn(2,1); end; for j=1:16; obsErrs(1:2,j) = randn(2,1); end; disp(evolErrs); disp(obsErrs);
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

double PI = 3.141592653589793;

double evolErrs_rowwise[] = {
-0.343003152130103,-0.766711794483284,-0.016814112314737, 0.684339759945504,-1.401783282955619,-1.521660304521858,-0.127785244107286, 0.602860572524585,-0.139677982915557, 0.407768714902350, 0.397539533883833,-0.317539749169638,-0.779285825610984,-1.935513755513929, 0.678730596165904,
1.666349045016822, 2.635481573310387, 0.304155468427342, 0.055808274805755,-1.360112379179931, 1.054743814037827,-1.410338023439304,-0.456929290517258,-0.983310072206319, 0.242994841538368,-0.175692485792199,-1.101615186229668,-1.762205119649466, 1.526915548584107,-2.277161011565906 };

double obsErrs_rowwise[] = {
-1.428567988496096, 0.913205695955837,-1.576872295738796,-1.888336147279610, 1.116853507009928, 1.615888145666843,-0.102585012191329,-0.192732954692481, 0.160906008337421,-0.024849020282298,-1.001561909251739,-0.314462113181954,0.276865687293751, 0.175430340572582, 0.746792737753047, 1.648965874319728,
-1.114618464565160, 0.976371425014641, 0.204080086636545, 0.736193913185726, 0.743379272133998,-1.666530392059792, 0.622727541956653, 0.794595441386172, 0.539084689771962,-2.548385761079745,-1.161623730001803, 1.066876935479899,1.748562141782206, 0.362976707912966, 0.842263598054067, 1.725578381396231
};
    

int main(int argc, char* argv[]) {

	printf("rotation starting\n");

	int i;
	
	double alpha = 2.0 * PI / 16.0;

	double F_rowwise[] = {
			cos(alpha), -sin(alpha),
			sin(alpha),  cos(alpha),
	};

	double G_rowwise[] = {
	  1, 0,
	  0, 1,
	  1, 1,
	  2, 1,
	  1, 2,
	  3, 1
	};

	double evolutionStd   = 1e-3;
	double observationStd = 1e-1;

	int k = 16;

	int obs_dim = 2;

	kalman_matrix_t* evolErrs = matrix_create_from_rowwise(evolErrs_rowwise, 2, 15);
	kalman_matrix_t* obsErrs  = matrix_create_from_rowwise(obsErrs_rowwise , 2, 16);

	kalman_matrix_t* F = matrix_create_from_rowwise(F_rowwise, 2, 2);
	kalman_matrix_t* G = matrix_create_from_rowwise(G_rowwise, 6, 2);
	
	printf("F = \n");
	matrix_print(F, "%.6f");

	G = matrix_create_sub(G,0,obs_dim,0,2);

	printf("G = \n");
	matrix_print(G, "%.6f");

	char K_type = 'w';
	kalman_matrix_t* K = matrix_create_constant(2,2,0.0);
	for (i=0; i<2; i++) matrix_set(K, i,i, evolutionStd);
	
	printf("K = \n");
	matrix_print(K, "%.3e");

	char C_type = 'w';
	kalman_matrix_t* C = matrix_create_constant(obs_dim,obs_dim,0.0);
	for (i=0; i<obs_dim; i++) matrix_set(C, i,i, observationStd);
	
	printf("C = \n");
	matrix_print(C, "%.3e");

	kalman_matrix_t* states = matrix_create_constant(2,      k, 0.0);
	kalman_matrix_t* obs    = matrix_create_constant(obs_dim,k, 0.0);

	matrix_set(states, 0, 0, 1.0);
	matrix_set(states, 1, 0, 0.0);

	for (i=1; i<k; i++) {
	  matrix_multiply_accumulate(states,   i, 0,
				     F,        0, 0,
				     states, i-1, 0,
				     2, 1, matrix_cols(F));

          matrix_multiply_accumulate(states,     i, 0,
				     K,          0, 0,
				     evolErrs, i-1, 0,
				     2, 1, matrix_cols(K));
	}

	printf("states = \n");
	matrix_print(states, "%.3f");

	for (i=0; i<k; i++) {
	  matrix_multiply_accumulate(obs,      i, 0,
				     G,        0, 0,
				     states,   i, 0,
				     2, 1, matrix_cols(G));

	  matrix_multiply_accumulate(states,     i, 0,
				     C,          0, 0,
				     obsErrs,    i, 0,
				     2, 1, matrix_cols(C));
	}

	printf("obs = \n");
	matrix_print(obs, "%.3f");
	
	//K = CovarianceMatrix(evolutionStd^-1  *ones(2,1),'w');
	//C = CovarianceMatrix(observationStd^-1*ones(obs_dim,1), 'w');


	return 0;
}

// C(i:i+m,j:j+n) += A * B(k:k+m,l:l+n)
void matrix_multiply_accumulate(
				kalman_matrix_t* C, int i, int j,
				kalman_matrix_t* A, int p, int q,
				kalman_matrix_t* B, int k, int l,
				int Csub_rows, int Csub_cols, int Asub_cols 
				) {
  int r,c,s;

  for (r=0; r<Csub_rows; r++) {
    for (c=0; c<Csub_cols; c++) {
      double x = matrix_get(C, i+r, j+c);
      for (s=0; s<Asub_cols; s++) {
        x += matrix_get(A,r+p,s+q) * matrix_get(B, k+s, j+l);
      }
      matrix_set(C, i+r, j+c, x);
    }
  }
}
     
