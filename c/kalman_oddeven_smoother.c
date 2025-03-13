/*
 * kalman_oddeven.c
 *
 * (C) Sivan Toledo and Shahaf Gargir, 2024-2025
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

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
#include "parallel.h"
#include "memory.h"
#include "assertions.h"

/******************************************************************************/
/* KALMAN STEPS                                                               */
/******************************************************************************/

typedef struct step_st {
	kalman_step_index_t step; // logical step number
	int32_t dimension;

	kalman_matrix_t* H;
	kalman_matrix_t* F;
	kalman_matrix_t* c;

	kalman_matrix_t* G;
	kalman_matrix_t* o;
	kalman_matrix_t* C;

	kalman_matrix_t* X;
	kalman_matrix_t* Y;
	kalman_matrix_t* Z;
	kalman_matrix_t* R;

	kalman_matrix_t* X_tilde;
	kalman_matrix_t* F_tilde;
	kalman_matrix_t* H_tilde;
	kalman_matrix_t* R_tilde;
	kalman_matrix_t* G_tilde;

	kalman_matrix_t* state;
	kalman_matrix_t* covariance;
} step_t;

static void* step_create() {
	step_t* s = malloc(sizeof(step_t));
	s->step      = -1;
	s->dimension = -1;


	s->H = NULL;
	s->F = NULL;
	s->c = NULL;

	s->G = NULL;
	s->o = NULL;
	s->C = NULL;

	s->X = NULL;
	s->Y = NULL;
	s->Z = NULL;
	s->R = NULL;

	s->X_tilde = NULL;
	s->F_tilde = NULL;
	s->H_tilde = NULL;
	s->R_tilde = NULL;
	s->G_tilde = NULL;

	s->state = NULL;
	s->covariance = NULL;

	assert( s != NULL );
	return s;
}

static void step_free(void* v) {
	step_t* s = (step_t*) v;
	
	matrix_free(s->H);
	matrix_free(s->F);
	matrix_free(s->c);

	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->X);
	matrix_free(s->Y);
	matrix_free(s->Z);
	matrix_free(s->R);

	matrix_free(s->X_tilde);
	matrix_free(s->F_tilde);
	matrix_free(s->H_tilde);
	matrix_free(s->R_tilde);
	matrix_free(s->G_tilde);

	matrix_free(s->state);
	matrix_free(s->covariance);

	free(s);
}

static void step_rollback(void* v) {
	step_t* s = (step_t*) v;
	
	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->X);
	matrix_free(s->Y);
	matrix_free(s->Z);
	matrix_free(s->R);

	matrix_free(s->X_tilde);
	matrix_free(s->F_tilde);
	matrix_free(s->H_tilde);
	matrix_free(s->R_tilde);
	matrix_free(s->G_tilde);

	matrix_free(s->state);
	matrix_free(s->covariance);
}

static kalman_step_index_t step_get_index(void* v) {
	return ((step_t*) v)->step;
}

static int32_t step_get_dimension(void* v) {
	return ((step_t*) v)->dimension;
}

static kalman_matrix_t* step_get_state(void* v) {
	return ((step_t*) v)->state;
}

static kalman_matrix_t* step_get_covariance(void* v) {
	return ((step_t*) v)->covariance;
}

static char step_get_covariance_type(void* v) {
	return 'C';
}

/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

/*
 * In this implementation, evolve and observe only weigh the inputs and store
 * them in the array of steps.
 * It might make more sense to just store and to weigh later (will allow
 * unification with associative evolve and observe!)
 */

static void evolve(kalman_t* kalman, int32_t n_i, matrix_t* H_i, matrix_t* F_i, matrix_t* c_i, matrix_t* K_i, char K_type) {
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
	//assert(K_i!=NULL);


	matrix_t* V_i_H_i = kalman_covariance_matrix_weigh(K_i,K_type,H_i);
	matrix_t* V_i_F_i = kalman_covariance_matrix_weigh(K_i,K_type,F_i);
	matrix_t* V_i_c_i = kalman_covariance_matrix_weigh(K_i,K_type,c_i);

	matrix_mutate_scale(V_i_F_i,-1.0);

	kalman_current->H	= V_i_H_i;
	kalman_current->F 	= V_i_F_i;
	kalman_current->c  = V_i_c_i;
}

static void observe(kalman_t* kalman, matrix_t* G_i, matrix_t* o_i, matrix_t* C_i, char C_type) {

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

	matrix_t* W_i_G_i = NULL;
	matrix_t* W_i_o_i = NULL;

	if (o_i != NULL) {
#ifdef BUILD_DEBUG_PRINTOUTS
		printf("C_i(%c) ",C_type);
		matrix_print(C_i,"%.3e");
		printf("G_i ");
		matrix_print(G_i,"%.3e");
		printf("o_i ");
		matrix_print(o_i,"%.3e");
#endif
	W_i_G_i = kalman_covariance_matrix_weigh(C_i,C_type,G_i);
	W_i_o_i = kalman_covariance_matrix_weigh(C_i,C_type,o_i);

	kalman_current->G  = W_i_G_i;
	kalman_current->o  = W_i_o_i;

#ifdef BUILD_DEBUG_PRINTOUTS
		printf("W_i_G_i ");
		matrix_print(W_i_G_i,"%.3e");
		printf("W_i_o_i ");
		matrix_print(W_i_o_i,"%.3e");
#endif
	}

	farray_append(kalman->steps,kalman->current);
	kalman->current = NULL;

#ifdef BUILD_DEBUG_PRINTOUTS
	printf("kalman_observe done\n");
#endif
}

/******************************************************************************/
/* ADDITIONAL MATRIX OPERATIONS                                               */
/******************************************************************************/

static void apply_Q_on_block_matrix (matrix_t* R, matrix_t* Q, matrix_t** upper, matrix_t** lower) {
	matrix_t* concat = matrix_create_vconcat(*upper, *lower);
	matrix_mutate_apply_qt(R, Q, concat);
	matrix_t *new_upper = matrix_create_sub(concat,0,matrix_rows(*upper), 0, matrix_cols(concat));
	matrix_t *new_lower = matrix_create_sub(concat,matrix_rows(*upper),matrix_rows(*lower), 0, matrix_cols(concat));

	matrix_free(*upper);
	matrix_free(*lower);
	matrix_free(concat);

	*upper = new_upper;
	*lower = new_lower;
}

static void free_and_assign(matrix_t** original, matrix_t* new){
	if (*original != NULL){
		matrix_free(*original);
	}
	*original = new;
}

/******************************************************************************/
/* FUNCIONS THAT CAN BE APPLIED IN PARALLEL TO AN ARRAY OF STEPS              */
/******************************************************************************/

//void assign_steps(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end) {
static void create_steps_array(void* kalman_v, void* steps_v, int length, size_t start, size_t end) {
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;
	
	for (size_t i = start; i < end; ++i) {
		steps[i] = farray_get(kalman->steps,i);
	}
}

//void G_F_to_R_tilde(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void G_F_to_R_tilde(void** steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (size_t j_ = start; j_ < end; ++j_) {
		size_t j = j_ * 2;

		step_t* step_i = steps[j];
		matrix_t* G_i = step_i->G;

		if (j == length - 1){
			matrix_t* o_i = step_i->o;

			matrix_t* R_tilde = matrix_create_copy(G_i);

			matrix_t* TAU = matrix_create_mutate_qr(R_tilde);
			matrix_mutate_apply_qt(R_tilde,TAU,o_i);

			step_i->R_tilde = R_tilde;
			matrix_mutate_triu(step_i->R_tilde);
			matrix_free(TAU);
		}else{

			step_t* step_ipo = steps[j + 1];
			matrix_t* o_i = step_i->o;

			matrix_t* c_ipo = step_ipo->c;

			matrix_t* F_ipo = step_ipo->F;
			matrix_t* H_ipo = step_ipo->H;

			matrix_t* concat = matrix_create_vconcat(G_i, F_ipo);

			matrix_t* Q, *R_tilde;

			R_tilde = concat;
			Q = matrix_create_mutate_qr(R_tilde);

			step_i->R_tilde = matrix_create_sub(R_tilde,0,matrix_cols(R_tilde), 0, matrix_cols(R_tilde));
			matrix_mutate_triu(step_i->R_tilde);
			step_i->X = matrix_create_constant(matrix_rows(G_i), matrix_cols(H_ipo), 0);

			free_and_assign(&(step_ipo->H_tilde), matrix_create_copy(H_ipo));

			apply_Q_on_block_matrix(R_tilde, Q, &(step_i->o), &(step_ipo->c));
			apply_Q_on_block_matrix(R_tilde, Q, &(step_i->X), &(step_ipo->H_tilde));
			matrix_free(R_tilde);
			matrix_free(Q);
		}
	}//Done first part of the algorithm
}

//void H_R_tilde_to_R(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void H_R_tilde_to_R(void** steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (size_t j_ = start; j_ < end; j_++){
		size_t j = j_ * 2;

		step_t* step_i = steps[j];
		if (j == 0){ //First index
			step_i->R = matrix_create_copy(step_i->R_tilde);
			continue;
		}

		matrix_t* R_tilde = step_i->R_tilde;
		matrix_t* F_i = step_i->F;
		matrix_t* H_i = step_i->H;

		matrix_t* concat = matrix_create_vconcat(H_i, R_tilde);
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));
		
		
		matrix_mutate_triu(R);
		step_i->R = R;

		step_i->Z = matrix_create_constant(matrix_rows(R_tilde), matrix_cols(F_i), 0);

		step_i->F_tilde = matrix_create_copy(F_i);
		apply_Q_on_block_matrix(R_tall, Q, &(step_i->F_tilde), &(step_i->Z));
		apply_Q_on_block_matrix(R_tall, Q, &(step_i->c), &(step_i->o));

		if(j + 1 != length){
			matrix_t *X = step_i->X;

			step_i->Y = matrix_create_constant(matrix_rows(F_i), matrix_cols(X), 0);

			step_i->X_tilde = matrix_create_copy(X);

			apply_Q_on_block_matrix(R_tall, Q, &(step_i->Y), &(step_i->X_tilde));

		}
		matrix_free(R_tall);
		matrix_free(Q);
	} //Done second part of the algorithm

}

//void H_tilde_G_to_G_tilde(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void H_tilde_G_to_G_tilde(void** steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (size_t j_ = start; j_ < end; j_++){
		size_t j = j_ * 2;
		// Sivan Feb 2025 wrong type and seens not to be used, commenting out
		//int i = steps[j];

		step_t* step_ipo = steps[j + 1];
		matrix_t* H_tilde = step_ipo->H_tilde;
		matrix_t* G_ipo = step_ipo->G;

		matrix_t* concat = matrix_create_vconcat(H_tilde, G_ipo);
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));

		matrix_mutate_triu(R);

		free_and_assign(&(step_ipo->G_tilde), R); 


		apply_Q_on_block_matrix(R_tall, Q, &(step_ipo->c), &(step_ipo->o));
		matrix_free(R_tall);
		matrix_free(Q);


	}//Done third part of the algorithm
}

//void Variables_Renaming(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void Variables_Renaming(void** steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (size_t j_ = start; j_ < end; ++j_){
		size_t j = j_ * 2;
		// int i = steps[j];
		// int ipo = steps[j + 1];
		step_t* step_ipo = steps[j + 1];
		matrix_t* G_tilde = step_ipo->G_tilde;
		matrix_t* o = step_ipo->c;

		matrix_t * copy_G_tilde = matrix_create_copy(G_tilde);
		matrix_t * copy_o = matrix_create_copy(o);

		if (j != 0){
			step_t* step_i = steps[j];
			matrix_t* Z = step_i->Z;
			matrix_t* X_tilde = step_i->X_tilde;
			matrix_t* c = step_i->o;

			matrix_t * copy_Z = matrix_create_copy(Z);
			matrix_t * copy_X_tilde = matrix_create_copy(X_tilde);
			matrix_t * copy_c = matrix_create_copy(c);

			free_and_assign(&(step_ipo->F), copy_Z); 

			free_and_assign(&(step_ipo->H), copy_X_tilde); 

			free_and_assign(&(step_ipo->c), copy_c);
		}
		free_and_assign(&(step_ipo->G), copy_G_tilde);
		free_and_assign(&(step_ipo->o), copy_o);
	}
}

//void Init_new_steps(void* new_steps_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void extract_recursion_steps(void* new_steps_v, void* steps_v, int length, size_t start, size_t end){
	step_t** new_steps = (step_t**) new_steps_v;
	step_t** steps     = (step_t**) steps_v;

    for (size_t i = start; i < end; ++i) {
		new_steps[i] = steps[2*i + 1];
	}
}

//void Solve_Estimates(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void Solve_Estimates(void* steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;

	for (size_t j_ = start; j_ < end; ++j_){
		size_t j = j_ * 2;

		step_t* step_i = steps[j];
		matrix_t* R = step_i->R;

		if (j == 0){

			step_t* step_ipo = steps[j + 1];
			matrix_t* x_ipo = step_ipo->state;

			matrix_t* X = step_i->X;
			matrix_t* o_i = step_i->o;
			matrix_t* mul = matrix_create_constant(matrix_rows(X), 1, 0);
			matrix_mutate_gemm(1, X, x_ipo, 0, mul);
			matrix_t* new_b = matrix_create_subtract(o_i, mul);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul);
			matrix_free(new_b);

		}else if (j != length - 1){
			
			step_t* step_ipo = steps[j + 1];
			matrix_t* x_ipo = step_ipo->state;
			step_t* step_imo = steps[j - 1];
			matrix_t* x_imo = step_imo->state;

			matrix_t* F_tilde = step_i->F_tilde;
			matrix_t* Y = step_i->Y;
			matrix_t* c = step_i->c;

			matrix_t* mul1 = matrix_create_constant(matrix_rows(F_tilde), 1, 0);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul1);
			matrix_t* mul2 = matrix_create_constant(matrix_rows(Y), 1, 0);
			matrix_mutate_gemm(1, Y, x_ipo, 0, mul2);
			matrix_t* new_b_mid = matrix_create_subtract(c, mul1);
			matrix_t* new_b = matrix_create_subtract(new_b_mid, mul2);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul1);
			matrix_free(mul2);
			matrix_free(new_b_mid);
			matrix_free(new_b);
		}else{

			
			step_t* step_imo = steps[j - 1];
			matrix_t* x_imo = step_imo->state;

			matrix_t* F_tilde = step_i->F_tilde;
			matrix_t* c = step_i->c;

			matrix_t* mul = matrix_create_constant(matrix_rows(F_tilde), 1, 0);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul);
			matrix_t* new_b = matrix_create_subtract(c, mul);

			matrix_mutate_triu(R);

			step_i->state = matrix_create_trisolve(R,new_b);

			matrix_free(mul);
			matrix_free(new_b);
		}
	}
}
// ==========================================
// Cov Change 4
// ==========================================

//void Convert_LDLT(void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void Convert_LDLT(void* steps_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;

	for (size_t j_ = start; j_ < end; ++j_){
		size_t j = j_ * 2;

		step_t* step = steps[j];
		matrix_t* R = step->R;
		matrix_t* R_inv = matrix_create_inverse(R);
		step->R = R_inv;
		
		if (j == 0){
			matrix_t* X = step->X;
			matrix_t* X_inv = matrix_create_trisolve(R,X);
			matrix_free(X);
			step->X = X_inv;

		}else if (j != length - 1){
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve(R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;	

			matrix_t* Y = step->Y;
			matrix_t* Y_inv = matrix_create_trisolve(R,Y);
			matrix_free(Y);
			step->Y = Y_inv;	
			
		}else{
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve(R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;		
		}
	}
}

//void SelInv(void* kalman_v, void* steps_v, int length, int** converters, size_t start, size_t end){
static void SelInv(void* steps_v, void** converters_v, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;
	int** converters = (int**) converters_v;

	int * vtop_n = converters[0];
	int * ptov_n2 = converters[1];
	
	for (size_t j_ = start; j_ < end; ++j_){
		size_t j = j_ * 2;

		step_t* step = steps[j];
	// 	Ainv(inz,j) = -Ainv(inz,inz) * L(inz,j);
	// 	Ainv(j,inz) = Ainv(inz,j)';
	// 	Ainv(j,j) = 1/D(j,j) - Ainv(j,inz) * L(inz,j);
		
		if (j == 0){
			int inz1 = 0;
			
			int physical_column = (int) ceil(((double)length)/2) + ptov_n2[inz1];
			matrix_t* L_j_inz = step->X;

			step_t *physical_step = steps[vtop_n[physical_column]];
			int physical_row = physical_step->step;

			matrix_t* Ainv_inz_inz = physical_step->R;

			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
			matrix_t* Dinv = step->R;

			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_create_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);
			matrix_free(Dinv_Dinv_t);

			matrix_free(step->X);
			step->X = Ainv_j_inz;

			matrix_free(step->R);
			step->R = Ainv_j_j;

		}else if (j != length - 1){
			int inz1 = j/2 - 1;
			int inz2 = j/2;

			int physical_column1 = ceil(((double)length)/2) + ptov_n2[inz1];

			int physical_column2 = ceil(((double)length)/2) + ptov_n2[inz2];

			int regular_order = 1;

	// 		We have two cases here: Y,F_tilde and F_tilde,Y 
	// 		and the order is matter.

			step_t *physical_step1;
			step_t *physical_step2;
			
			matrix_t* L_j_inz;

			if (physical_column1 < physical_column2){ //This is the ragular order

				matrix_t* F_tilde_t = matrix_create_transpose(step->F_tilde);
				matrix_t* Y_t = matrix_create_transpose(step->Y);
				matrix_t* L_j_inz_t = matrix_create_vconcat(F_tilde_t, Y_t);
				L_j_inz = matrix_create_transpose(L_j_inz_t);

				matrix_free(F_tilde_t);
				matrix_free(Y_t);
				matrix_free(L_j_inz_t);

				physical_step1 = steps[vtop_n[physical_column1]];
				int physical_row1 = physical_step1->step;
				physical_step2 = steps[vtop_n[physical_column2]];
				int physical_row2 = physical_step2->step;
			}else{
				regular_order = 0;

				matrix_t* F_tilde_t = matrix_create_transpose(step->F_tilde);
				matrix_t* Y_t = matrix_create_transpose(step->Y);
				matrix_t* L_j_inz_t = matrix_create_vconcat(Y_t, F_tilde_t);
				L_j_inz = matrix_create_transpose(L_j_inz_t);

				matrix_free(F_tilde_t);
				matrix_free(Y_t);
				matrix_free(L_j_inz_t);
					
				int tmp = physical_column1;
				physical_column1 = physical_column2;
				physical_column2 = tmp;

				physical_step1 = steps[vtop_n[physical_column1]];
				int physical_row1 = physical_step1->step;
				physical_step2 = steps[vtop_n[physical_column2]];
				int physical_row2 = physical_step2->step;	
			}

			
	// 		Now we have to do the same process in row 2 but
	// 		column 2. but we can't call virtual_physical

			int start_index = ceil(((double)length)/2);

			matrix_t* top_right;
	// 		again we have 3 cases:
			if (start_index == physical_column1){
				top_right = physical_step1->X;
			}else if (start_index + ceil(floor(((double)length)/2)/2) > physical_column1 + 1 || (length/2)%2 == 0){
	
				int source_index = floor(((double)vtop_n[physical_column2])/4);

	// 			Now lets check the source of F_tilde and Y of this row
				int F_tilde_location = physical_column1 - start_index - 1;
				if (source_index == F_tilde_location){	
					top_right = physical_step1->F_tilde;
				}else{
					top_right = physical_step1->Y;
				    assert(source_index == physical_column1 - start_index);
				}
			}else{	
				top_right = physical_step1->F_tilde;
			}



	// 		Ainv_inz_inz = [kalman.steps{physical_row1 + 1}.R , top_right
	// 						top_right'                        , kalman.steps{physical_row2 + 1}.R];
			matrix_t* top_right_t = matrix_create_transpose(top_right);
			matrix_t* R_top_right_t = matrix_create_vconcat(physical_step1->R, top_right_t);
			matrix_t* top_right_R = matrix_create_vconcat(top_right, physical_step2->R);
			matrix_t* t1 = matrix_create_transpose(R_top_right_t);
			matrix_t* t2 = matrix_create_transpose(top_right_R);
			matrix_t* Ainv_inz_inz_t = matrix_create_vconcat(t1, t2);
			matrix_t* Ainv_inz_inz = matrix_create_transpose(Ainv_inz_inz_t);

			matrix_free(top_right_t);
			matrix_free(R_top_right_t);
			matrix_free(top_right_R);
			matrix_free(t1);
			matrix_free(t2);
			matrix_free(Ainv_inz_inz_t);

	// 		Ainv_j_inz = -L_j_inz * Ainv_inz_inz;
			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_create_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);
			
			if (regular_order){ //This is the ragular order

				int F_tilde_cols = matrix_cols(step->F_tilde);
				int Y_cols = matrix_cols(step->Y);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, F_tilde_cols);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), F_tilde_cols, matrix_cols(Ainv_j_inz) - F_tilde_cols);	
			}else{
				int F_tilde_cols = matrix_cols(step->F_tilde);
				int Y_cols = matrix_cols(step->Y);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, Y_cols);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), Y_cols, matrix_cols(Ainv_j_inz) - Y_cols);
			}
			
			matrix_free(step->R);
			step->R = Ainv_j_j;

		}else{
			int inz1 = j/2 - 1;			

			int physical_column = ceil(((double)length)/2) + ptov_n2[inz1];

			matrix_t* L_j_inz = step->F_tilde;

			step_t *physical_step = steps[vtop_n[physical_column]];
			int physical_row = physical_step->step;

			matrix_t* Ainv_inz_inz = physical_step->R;

	// 		Ainv_j_inz = -L_j_inz * Ainv_inz_inz;
			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_create_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);

	// 		kalman.steps{i}.F_tilde = Ainv_j_inz;
			matrix_free(step->F_tilde);
			step->F_tilde = Ainv_j_inz;
	// 		kalman.steps{i}.R = Ainv_j_j;
			matrix_free(step->R);
			step->R = Ainv_j_j;
		}
		step->covariance = matrix_create_copy(step->R);
	}
}

// ==========================================
// Cov End Change 4
// ==========================================

// ==========================================
// Cov Change 2
// ==========================================


//void one_layer_converter (void* kalman_v, void* steps_v, int length, int** helper, size_t start, size_t end){
static void one_layer_converter (void* array, int length, size_t start, size_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	//step_t* *steps = (step_t**)steps_v;
	
	int** helper = (int**) array;
	 
	int* virtual_to_physical = helper[0];
	int* physical_to_virtual = helper[1];
	int* variables = helper[2];

	int n = variables[0];
	int start_ = variables[1];
	int jump = variables[2];
	int skip = variables[3];
	int index = variables[4];

	for (int j_ = start; j_ < end; ++j_){
		int i = j_ * jump + start_;

		virtual_to_physical[index + j_] = i;
		if (!skip){
			physical_to_virtual[i/2] = (index + j_) - ceil(((double)n)/2);
		}
	}
	variables[4] = index;
}

static int** virtual_physical(int n) {

	int* virtual_to_physical = (int*) malloc(    n * sizeof(int));
	int* physical_to_virtual = (int*) malloc((n/2) * sizeof(int));
	
	int* variables = (int *)malloc(5 * sizeof(int));

	int start = 0;
	int jump = 2;

	int index = 0;

	int skip = 1;

	int ** result = (int **)malloc(3 * sizeof(int *));
	result[0] = virtual_to_physical;
	result[1] = physical_to_virtual;
	result[2] = variables;

	while (start < n){

		variables[0] = n;
		variables[1] = start;
		variables[2] = jump;
		variables[3] = skip;
		variables[4] = index;

		//#ifdef PARALLEL
		//parallel_for_c(NULL, NULL, 0, result, ceil((double)(n - start)/jump),BLOCKSIZE,one_layer_converter);
		//#else
		//one_layer_converter(NULL, NULL, 0, result, 0, ceil((double)(n - start)/jump));
		//#endif
		foreach_in_range(one_layer_converter, result, ceil((double)(n - start)/jump) /* ??? */, ceil((double)(n - start)/jump));

		index += ceil((double)(n - start)/jump);

		skip = 0;

		start = start * 2 + 1;
		jump = jump * 2;
	}
	
	free(variables);
	return result;
}
// ==========================================
// End Cov Change 2
// ==========================================

//void smooth_recursive(kalman_t* kalman, step_t* *steps, int length) {
static void smooth_recursive(kalman_options_t options, step_t** steps, int length) {
	
    if (length == 1) {

		step_t* singleStep = steps[0];

		matrix_t* G = matrix_create_copy(singleStep->G);
		matrix_t* o = matrix_create_copy(singleStep->o);

		matrix_t* TAU = matrix_create_mutate_qr(G);
		matrix_mutate_apply_qt(G,TAU,o);

		matrix_mutate_triu(G);

		// ==========================================
		// Cov Change 3
		// ==========================================

		matrix_t* RT = matrix_create_transpose(G);

		matrix_t* RT_R = matrix_create_constant(matrix_cols(G), matrix_cols(G), 0);
		matrix_mutate_gemm(1, RT, G, 0, RT_R);

		singleStep->R = matrix_create_inverse(RT_R);

		singleStep->covariance = matrix_create_copy(singleStep->R);

		matrix_free(RT);
		matrix_free(RT_R);
		
		// ==========================================
		// End Cov Change 3
		// ==========================================
		singleStep->state = matrix_create_trisolve(G,o);
		
		matrix_free(TAU);
		matrix_free(G);
		matrix_free(o);
		
		return;
    }
	// First part of the algorithm
	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, (length + 1)/2,BLOCKSIZE,G_F_to_R_tilde);
	//#else
	//G_F_to_R_tilde(NULL, steps, length, NULL, 0, (length + 1)/2);
	//#endif

	foreach_in_range(G_F_to_R_tilde, steps, length, (length + 1)/2);

	//Second part of the algorithm
	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, (length + 1)/2,BLOCKSIZE, H_R_tilde_to_R);
	//#else
	//H_R_tilde_to_R(NULL, steps, length, NULL, 0, (length + 1)/2);
	//#endif
	foreach_in_range(H_R_tilde_to_R, steps, length, (length + 1)/2);

	//Third part of the algorithm
	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, length/2,BLOCKSIZE, H_tilde_G_to_G_tilde);
	//#else
	//H_tilde_G_to_G_tilde(NULL, steps, length, NULL, 0, length/2);
	//#endif
	// Sivan March 2025 not sure why this goes to length/2, not as in previous two
	foreach_in_range(H_tilde_G_to_G_tilde, steps, length, length/2);
	
	//Fix the last index

	if (length % 2 == 1){
		step_t* step_lmo = steps[length - 2];
		matrix_t* G_tilde = step_lmo->G_tilde;
		matrix_t* o_i = step_lmo->c;

		step_t* step_last = steps[length - 1];
		matrix_t* Z = step_last->Z;
		
		matrix_t* c_i = step_last->o;
		matrix_t* concat = matrix_create_vconcat(G_tilde, Z);
		
		matrix_t* R_tall = concat;
		matrix_t* Q = matrix_create_mutate_qr(R_tall);
		matrix_t* R = matrix_create_sub(R_tall,0,matrix_cols(R_tall), 0, matrix_cols(R_tall));

		matrix_mutate_triu(R);
		free_and_assign(&(step_lmo->G_tilde), R); //step_lmo->G_tilde = R;

		apply_Q_on_block_matrix(R_tall, Q, &(step_lmo->c), &(step_last->o));
		matrix_free(R_tall);
		matrix_free(Q);
	}

	// Create new Kalman for the recursion
	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, length/2, BLOCKSIZE, Variables_Renaming);
	//#else
	//Variables_Renaming(NULL, steps, length, NULL, 0, length/2);
	//#endif
	foreach_in_range(Variables_Renaming, steps, length, length/2);
	
	// The Recursion

	step_t** recursion_steps = (step_t**) malloc(length/2 * sizeof(step_t*));

	//#ifdef PARALLEL
	//parallel_for_c(new_steps, steps, length, NULL, length/2,BLOCKSIZE, Init_new_steps);
	//#else
	//Init_new_steps(new_steps, steps, length, NULL, 0, length/2);
	//#endif
	foreach_in_range_two(extract_recursion_steps, recursion_steps, steps, length, length/2);
	
	
	//smooth_recursive(NULL, new_steps, length/2);
	smooth_recursive(options, recursion_steps, length/2);

	free(recursion_steps);

	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, (length + 1)/2, BLOCKSIZE, Solve_Estimates);
	//#else
	//Solve_Estimates(NULL, steps, length, NULL, 0, (length + 1)/2);
	//#endif
	foreach_in_range(Solve_Estimates, steps, length, (length + 1)/2);

	// ==========================================
	// Change 1
	// ==========================================

	if ((options & KALMAN_NO_COVARIANCE) == 0) {
	// Since the selinv algorithm works with LDL^T matrices, we
	// start with converting our RR^T to LDL^T

	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, NULL, (length + 1)/2, BLOCKSIZE, Convert_LDLT);
	//#else
	//Convert_LDLT(NULL, steps, length, NULL, 0, (length + 1)/2);
	//#endif
	foreach_in_range(Convert_LDLT, steps, length, (length + 1)/2);

	// Now we can start the selinv algrithm for out case
	
	int** result = virtual_physical(length);

	//#ifdef PARALLEL
	//parallel_for_c(NULL, steps, length, result, (length + 1)/2, BLOCKSIZE, SelInv);
	//#else
	//SelInv(NULL, steps, length, result, 0, (length + 1)/2);
	//#endif
	foreach_in_range_two(SelInv, steps, result, length, (length + 1)/2);
	
	free(result[0]);
	free(result[1]);
	free(result);
	
	}
	// % ==========================================
	// % End Change 1
	// % ==========================================
	

}

static void smooth(kalman_t* kalman) {

	int length = farray_size(kalman->steps);
	step_t** steps = (step_t**) malloc(length * sizeof(step_t*));
	
//#ifdef PARALLEL
//	parallel_for_c(kalman, steps, length, NULL, length, BLOCKSIZE, assign_steps);
//#else
//	assign_steps(kalman, steps, length, NULL, 0, length);
//#endif
	foreach_in_range_two(create_steps_array, kalman, steps, length, length);

	//for (int i = 0; i < length; ++i) {
	//	steps[i] = farray_get(kalman->steps,i);
	//}

	//parallel_smooth(kalman, steps, length);
	//smooth_recursive(NULL, steps, length);
	smooth_recursive(kalman->options, steps, length);

	free(steps);
}

void kalman_create_oddeven(kalman_t* kalman) {
	kalman->evolve  = evolve;	
	kalman->observe = observe;	
	kalman->smooth  = smooth;	
	
	kalman->step_create        = step_create;
	kalman->step_free          = step_free;
	kalman->step_rollback      = step_rollback;
	kalman->step_get_index     = step_get_index;
	kalman->step_get_dimension = step_get_dimension;
	kalman->step_get_state     = step_get_state;
	kalman->step_get_covariance      = step_get_covariance;
	kalman->step_get_covariance_type = step_get_covariance_type;
}
/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
