/*
 * kalman_oddeven.c
 *
 * A parallel linear Kalman smoother as propsed in the article:
 *  Shahaf Gargir and Sivan Toledo.
 * Parallel-in-time Kalman smoothing using orthogonal transformations.
 * In Proceedings of the IEEE International Parallel and Distributed Processing Symposium (IPDPS), 2025.
 * https://arxiv.org/abs/2502.11686
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
#include "memory.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

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

#ifdef OBSOLETE
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
#endif
/******************************************************************************/
/* KALMAN                                                                     */
/******************************************************************************/

/*
 * In this implementation, evolve and observe only weigh the inputs and store
 * them in the array of steps.
 * It might make more sense to just store and to weigh later (will allow
 * unification with associative evolve and observe!)
 */

#ifdef OBSOLETE
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
#endif

/******************************************************************************/
/* ADDITIONAL MATRIX OPERATIONS                                               */
/******************************************************************************/

static void apply_Q_on_block_matrix (matrix_t* R, matrix_t* Q, matrix_t** upper, matrix_t** lower) {
  //printf(">>> %d\n",__LINE__);
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

/*
 * Apply QT to a 2-by-1 block matrix, splitting the transformed
 * concatenation at rows_top.
 */
static void apply_QT_to_block_matrix (matrix_t* QR, matrix_t* TAU, int32_t rows_top, matrix_t** upper, matrix_t** lower) {
  //printf(">>> %d\n",__LINE__);
    matrix_t* concat = matrix_create_vconcat(*upper, *lower);
    matrix_mutate_apply_qt(QR, TAU, concat);
    //printf("rows_top=%d rest=%d total=%d\n",rows_top,matrix_rows(concat)-rows_top,matrix_rows(concat));
    matrix_t *new_upper = matrix_create_sub(concat,0,rows_top, 0, matrix_cols(concat));
    matrix_t *new_lower = matrix_create_sub(concat,rows_top,matrix_rows(concat)-rows_top, 0, matrix_cols(concat));

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

static void create_steps_array(void* kalman_v, void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end) {
    kalman_t* kalman = (kalman_t*) kalman_v;
    step_t** steps = (step_t**) steps_v;

    for (kalman_step_index_t i = start; i < end; ++i) {
        steps[i] = farray_get(kalman->steps,i);
    }
}

//void assign_steps(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end) {
static void create_steps_array_new(void* kalman_v, void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end) {
	kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;
	
	for (kalman_step_index_t i = start; i < end; ++i) {
		steps[i] = farray_get(kalman->steps,i);
	}
}

//void G_F_to_R_tilde(void* kalman_v, void* steps_v, kalman_step_index_t length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void G_F_to_R_tilde(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; ++j_) {
		kalman_step_index_t j = j_ * 2;

		step_t* step_i = steps[j];
		matrix_t* G_i = step_i->G;

		if (j == length - 1) {
			matrix_t* o_i = step_i->o;

            //printf(">>> %d\n",__LINE__);
			matrix_t* QR = matrix_create_copy(G_i);

	          //printf(">>> %d\n",__LINE__);
			matrix_t* TAU = matrix_create_mutate_qr(QR);
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_apply_qt(QR,TAU,o_i);

			step_i->R_tilde = QR;
			matrix_mutate_triu(step_i->R_tilde);
            //printf("~~~ A set R_tilde step %d dim %d R_tilde %d-by-%d\n",step_i->step,step_i->dimension,matrix_rows(step_i->R_tilde),matrix_cols(step_i->R_tilde));

			matrix_free(TAU);
		} else {

			step_t* step_ipo = steps[j + 1];
			matrix_t* o_i = step_i->o;

			matrix_t* c_ipo = step_ipo->c;

			matrix_t* F_ipo = step_ipo->F;
			matrix_t* H_ipo = step_ipo->H;

            //printf(">>> %d\n",__LINE__);
            //printf("~~~ D set R_tilde step %d dim %d Gi? %d Fipo? %d Hipo? %d\n",step_i->step,step_i->dimension,G_i!=NULL,F_ipo!=NULL,H_ipo!=NULL);
			matrix_t* QR = matrix_create_vconcat(G_i, F_ipo);
            //printf(">>> %d\n",__LINE__);

            int nn,mm,ll;

            mm = matrix_rows(G_i);
            nn = matrix_cols(F_ipo);
            //printf(">>> %d\n",__LINE__);

			matrix_t* TAU;
			TAU = matrix_create_mutate_qr(QR);

            ll = MIN(nn,matrix_rows(QR));
            //printf("~~~ C set R_tilde step %d dim %d rows(QR)=%d rows(G_i)=%d rows(H_ipo)=%d\n",step_i->step,step_i->dimension,matrix_rows(QR),matrix_rows(G_i),matrix_rows(H_ipo));

            //step_i->R_tilde = matrix_create_sub(R_tilde,0,matrix_cols(R_tilde), 0, matrix_cols(R_tilde));
            step_i->R_tilde = matrix_create_sub(QR,0,ll, 0, nn);
			matrix_mutate_triu(step_i->R_tilde);
            //printf("~~~ B set R_tilde step %d dim %d R_tilde %d-by-%d\n",step_i->step,step_i->dimension,matrix_rows(step_i->R_tilde),matrix_cols(step_i->R_tilde));

            //step_i->X = matrix_create_constant(matrix_rows(G_i), matrix_cols(H_ipo), 0);
            // seems to cause a segmentation fault!!!
            step_i->X = matrix_create_constant(matrix_rows(QR)-matrix_rows(H_ipo), matrix_cols(H_ipo), 0); // number of rows is for the current split, above H_ipo
            //step_i->X = matrix_create_constant(matrix_rows(QR)-ll, matrix_cols(H_ipo), 0); // number of rows is for the current split, above H_ipo

			free_and_assign(&(step_ipo->H_tilde), matrix_create_copy(H_ipo));
            //printf(",,, step %d rows(H_tilde)=%d rows(X)=%d rows(c)=%d (added %d zero rows, should be %d)\n",step_ipo->step,matrix_rows(step_ipo->H_tilde),matrix_rows(step_i->X),matrix_rows(step_ipo->c),matrix_rows(QR)-ll,matrix_rows(QR)-matrix_rows(H_ipo));

			  //printf(">>> %d\n",__LINE__);
			  //printf("size(G_i)=%d,%d size(F_ipo)=%d,%d step=%d\n",matrix_rows(G_i),matrix_cols(G_i),matrix_rows(F_ipo),matrix_cols(F_ipo),step_ipo->step);
              //printf("mm=%d nn=%d rows(QR)=%d, cols(QR)=%d\n",mm,nn,matrix_rows(QR),matrix_cols(QR));
              //printf("ll=%d len(o)=%d len(c)=%d, rows(X)=%d rows(Htilde)=%d rows(QR)=%d\n",ll,matrix_rows(step_i->o),matrix_rows(step_ipo->c),matrix_rows(step_i->X),matrix_rows(step_ipo->H_tilde),matrix_rows(QR));
			apply_QT_to_block_matrix(QR, TAU, ll, &(step_i->o), &(step_ipo->c));
			  //printf(">>> %d setting X, o for step %d\n",__LINE__,step_i->step);
			apply_QT_to_block_matrix(QR, TAU, ll, &(step_i->X), &(step_ipo->H_tilde));
			//printf("... step %d rows(H_tilde)=%d rows(X)=%d rows(c)=%d\n",step_ipo->step,matrix_rows(step_ipo->H_tilde),matrix_rows(step_i->X),matrix_rows(step_ipo->c));
			matrix_free(QR);
			matrix_free(TAU);
		}
	}//Done first part of the algorithm
}

//void H_R_tilde_to_R(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void H_R_tilde_to_R(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; j_++){
		kalman_step_index_t j = j_ * 2;

		step_t* step_i = steps[j];
		if (j == 0){ //First index
			step_i->R = matrix_create_copy(step_i->R_tilde);
	        //printf("--- A set R step %d dim %d R %d-by-%d\n",step_i->step,step_i->dimension,matrix_rows(step_i->R),matrix_cols(step_i->R));

			continue;
		}

		matrix_t* R_tilde = step_i->R_tilde;
		matrix_t* F_i = step_i->F;
		matrix_t* H_i = step_i->H;

		int32_t nn,ll;
        ll = matrix_rows(H_i);
		nn = matrix_cols(H_i);
        //printf(">>> %d\n",__LINE__);

		matrix_t* QR  = matrix_create_vconcat(H_i, R_tilde);
		matrix_t* TAU = matrix_create_mutate_qr(QR);

		matrix_t* R   = matrix_create_sub(QR,0,nn, 0, nn);
		matrix_mutate_triu(R);
		step_i->R = R;
        //printf("--- B set R step %d dim %d R %d-by-%d\n",step_i->step,step_i->dimension,matrix_rows(step_i->R),matrix_cols(step_i->R));
        //printf("--- nn=%d QR %d-by=%d\n",nn,matrix_rows(QR),matrix_cols(QR));


		step_i->Z = matrix_create_constant(matrix_rows(R_tilde), matrix_cols(F_i), 0);
		step_i->F_tilde = matrix_create_copy(F_i);

		// the split should be after MIN(nn,matrix_rows(QR)), in case R is rectangular; in Matlab this is still nn!
		int32_t mm = MIN(nn, matrix_rows(QR));
		  //printf(">>> %d\n",__LINE__);
          //printf("   top=%d bottom=%d cols=%d split=%d oldtop=%d zerofill=%d\n",nn,matrix_rows(R_tilde),ll,mm,matrix_rows(step_i->F_tilde),matrix_rows(step_i->Z));
          //printf("   top v=%d bottom v=%d\n",matrix_rows(step_i->c),matrix_rows(step_i->o));
		apply_QT_to_block_matrix(QR, TAU, mm, &(step_i->F_tilde), &(step_i->Z));
		  //printf(">>> %d\n",__LINE__);
		apply_QT_to_block_matrix(QR, TAU, mm, &(step_i->c), &(step_i->o));
        //printf(">>> %d\n",__LINE__);

		if(j + 1 != length){
			matrix_t *X = step_i->X;

			step_i->Y = matrix_create_constant(matrix_rows(QR)-matrix_rows(X), matrix_cols(X), 0);

			step_i->X_tilde = matrix_create_copy(X);

			  //printf(">>> %d setting X for step %d (but not o)\n",__LINE__,step_i->step);
	          //printf("   top=%d bottom=%d split=%d\n",matrix_rows(step_i->Y),matrix_rows(step_i->X_tilde),mm);
			apply_QT_to_block_matrix(QR, TAU, mm, &(step_i->Y), &(step_i->X_tilde));

		}
		matrix_free(QR);
		matrix_free(TAU);
	} //Done second part of the algorithm

}

//void H_tilde_G_to_G_tilde(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void H_tilde_G_to_G_tilde(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; j_++){
		kalman_step_index_t j = j_ * 2;
		// Sivan Feb 2025 wrong type and seens not to be used, commenting out
		//int i = steps[j];

		step_t* step_ipo = steps[j + 1];
		matrix_t* H_tilde = step_ipo->H_tilde;
		matrix_t* G_ipo = step_ipo->G;

        //printf(">>> %d step %d rows(H_tilde)=%d rows(c)=%d\n",__LINE__,step_ipo->step,matrix_rows(H_tilde),matrix_rows(step_ipo->c));

        //printf("   top=%d bottom=%d cols=%d\n",matrix_rows(H_tilde),matrix_rows(G_ipo),matrix_cols(H_tilde));
		matrix_t* QR  = matrix_create_vconcat(H_tilde, G_ipo);
		matrix_t* TAU = matrix_create_mutate_qr(QR);

		int32_t ll = MIN(matrix_rows(QR),matrix_cols(QR));
		matrix_t* R   = matrix_create_sub(QR,0,ll, 0, matrix_cols(QR)); // R might be rectangular
		matrix_mutate_triu(R);

		free_and_assign(&(step_ipo->G_tilde), R); 

        //int32_t nn = MIN(ll,matrix_cols(H_tilde));
		  //printf(">>> %d\n",__LINE__);
	      //printf("   top=%d bottom=%d split=%d QR %d-by-%d\n",matrix_rows(step_ipo->c),matrix_rows(step_ipo->o),ll,matrix_rows(QR),matrix_cols(QR));
		apply_QT_to_block_matrix(QR, TAU, ll, &(step_ipo->c), &(step_ipo->o));
        //printf(">>> %d\n",__LINE__);
		matrix_free(QR);
		matrix_free(TAU);

	}//Done third part of the algorithm
}

//void Variables_Renaming(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void Variables_Renaming(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t* *steps = (step_t**)steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; ++j_){
		kalman_step_index_t j = j_ * 2;
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

//void Init_new_steps(void* new_steps_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void extract_recursion_steps(void* new_steps_v, void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	step_t** new_steps = (step_t**) new_steps_v;
	step_t** steps     = (step_t**) steps_v;

    for (kalman_step_index_t i = start; i < end; ++i) {
		new_steps[i] = steps[2*i + 1];
	}
}

//void Solve_Estimates(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void Solve_Estimates(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; ++j_){
		kalman_step_index_t j = j_ * 2;

		step_t* step_i = steps[j];
		matrix_t* R = step_i->R;

		if (j == 0){

			step_t* step_ipo = steps[j + 1];
			matrix_t* x_ipo = step_ipo->state;

			matrix_t* X = step_i->X;
			matrix_t* o_i = step_i->o;
			matrix_t* mul = matrix_create_constant(matrix_rows(X), 1, 0);
			  //printf(">>> %d R0 X0\n",__LINE__);
			  //printf("   step=%d\n",step_i->step);
              //printf("   gemm X_i %d-by-%d x_ipo(%d) %d-by-%d -> %d-by=%d\n",matrix_rows(X),matrix_cols(X),step_ipo->step,matrix_rows(x_ipo),matrix_cols(x_ipo),matrix_rows(mul),matrix_cols(mul));
              //printf("   subtract from o_i %d-by-%d solve with R %d-by-%d\n",matrix_rows(o_i),matrix_cols(o_i),matrix_rows(R),matrix_cols(R));
			matrix_mutate_gemm(1, X, x_ipo, 0, mul);
			matrix_t* new_b = matrix_create_subtract(o_i, mul);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve("U",R,new_b);
			//printf("   1 set state %d actual %d dim %d\n",step_i->step,matrix_rows(step_i->state),step_i->dimension);

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
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul1);
			matrix_t* mul2 = matrix_create_constant(matrix_rows(Y), 1, 0);
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, Y, x_ipo, 0, mul2);
			matrix_t* new_b_mid = matrix_create_subtract(c, mul1);
			matrix_t* new_b = matrix_create_subtract(new_b_mid, mul2);

			matrix_mutate_triu(R);
			step_i->state = matrix_create_trisolve("U",R,new_b);
            //printf("   2 set state %d actual %d dim %d\n",step_i->step,matrix_rows(step_i->state),step_i->dimension);

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
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, F_tilde, x_imo, 0, mul);
			matrix_t* new_b = matrix_create_subtract(c, mul);

			matrix_mutate_triu(R);

			step_i->state = matrix_create_trisolve("U",R,new_b);
            //printf("   3 set state %d actual %d dim %d\n",step_i->step,matrix_rows(step_i->state),step_i->dimension);

			matrix_free(mul);
			matrix_free(new_b);
		}
	}
}
// ==========================================
// Cov Change 4
// ==========================================

//void Convert_LDLT(void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void Convert_LDLT(void* steps_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;

	for (kalman_step_index_t j_ = start; j_ < end; ++j_){
		kalman_step_index_t j = j_ * 2;

		step_t* step = steps[j];
		matrix_t* R = step->R;
		matrix_t* R_inv = matrix_create_inverse(R);
		step->R = R_inv;
		
		if (j == 0){
			matrix_t* X = step->X;
			matrix_t* X_inv = matrix_create_trisolve("U",R,X);
			matrix_free(X);
			step->X = X_inv;

		}else if (j != length - 1){
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve("U",R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;	

			matrix_t* Y = step->Y;
			matrix_t* Y_inv = matrix_create_trisolve("U",R,Y);
			matrix_free(Y);
			step->Y = Y_inv;	
			
		}else{
			matrix_t* F_tilde = step->F_tilde;
			matrix_t* F_tilde_inv = matrix_create_trisolve("U",R,F_tilde);
			matrix_free(F_tilde);
			step->F_tilde = F_tilde_inv;		
		}
	}
}

//void SelInv(void* kalman_v, void* steps_v, int length, int** converters, kalman_step_index_t start, kalman_step_index_t end){
static void SelInv(void* steps_v, void* converters_v, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	step_t** steps = (step_t**) steps_v;
	kalman_step_index_t** converters = (kalman_step_index_t**) converters_v;

	kalman_step_index_t * vtop_n = converters[0];
	kalman_step_index_t * ptov_n2 = converters[1];
	
	for (kalman_step_index_t j_ = start; j_ < end; ++j_){
		kalman_step_index_t j = j_ * 2;

		step_t* step = steps[j];
	// 	Ainv(inz,j) = -Ainv(inz,inz) * L(inz,j);
	// 	Ainv(j,inz) = Ainv(inz,j)';
	// 	Ainv(j,j) = 1/D(j,j) - Ainv(j,inz) * L(inz,j);
		
		if (j == 0){
			kalman_step_index_t inz1 = 0;
			
			kalman_step_index_t physical_column = (kalman_step_index_t) ceil(((double)length)/2) + ptov_n2[inz1];
			matrix_t* L_j_inz = step->X;

			step_t *physical_step = steps[vtop_n[physical_column]];
			kalman_step_index_t physical_row = physical_step->step;

			matrix_t* Ainv_inz_inz = physical_step->R;

			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
			matrix_t* Dinv = step->R;

			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			  //printf(">>> %d\n",__LINE__);
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
			kalman_step_index_t inz1 = j/2 - 1;
			kalman_step_index_t inz2 = j/2;

			kalman_step_index_t physical_column1 = (kalman_step_index_t) ceil(((double)length)/2) + ptov_n2[inz1];

			kalman_step_index_t physical_column2 = (kalman_step_index_t) ceil(((double)length)/2) + ptov_n2[inz2];

			kalman_step_index_t regular_order = 1;

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
				kalman_step_index_t physical_row1 = physical_step1->step;
				physical_step2 = steps[vtop_n[physical_column2]];
				kalman_step_index_t physical_row2 = physical_step2->step;
			}else{
				regular_order = 0;

				matrix_t* F_tilde_t = matrix_create_transpose(step->F_tilde);
				matrix_t* Y_t = matrix_create_transpose(step->Y);
				matrix_t* L_j_inz_t = matrix_create_vconcat(Y_t, F_tilde_t);
				L_j_inz = matrix_create_transpose(L_j_inz_t);

				matrix_free(F_tilde_t);
				matrix_free(Y_t);
				matrix_free(L_j_inz_t);
					
				kalman_step_index_t tmp = physical_column1;
				physical_column1 = physical_column2;
				physical_column2 = tmp;

				physical_step1 = steps[vtop_n[physical_column1]];
				kalman_step_index_t physical_row1 = physical_step1->step;
				physical_step2 = steps[vtop_n[physical_column2]];
				kalman_step_index_t physical_row2 = physical_step2->step;
			}

			
	// 		Now we have to do the same process in row 2 but
	// 		column 2. but we can't call virtual_physical

			kalman_step_index_t start_index = (kalman_step_index_t) ceil(((double)length)/2);

			matrix_t* top_right;
	// 		again we have 3 cases:
			if (start_index == physical_column1){
				top_right = physical_step1->X;
			}else if (start_index + ceil(floor(((double)length)/2)/2) > physical_column1 + 1 || (length/2)%2 == 0){
	
				kalman_step_index_t source_index = (kalman_step_index_t) floor(((double)vtop_n[physical_column2])/4);

	// 			Now lets check the source of F_tilde and Y of this row
				kalman_step_index_t F_tilde_location = physical_column1 - start_index - 1;
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
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, Ainv_j_inz, L_j_inz_t, 0, tmp);

			matrix_t* Ainv_j_j = matrix_create_subtract(Dinv_Dinv_t, tmp);
			matrix_free(tmp);
			matrix_free(L_j_inz_t);
			
			if (regular_order){ //This is the ragular order

				int32_t F_tilde_cols = matrix_cols(step->F_tilde);
				int32_t Y_cols = matrix_cols(step->Y);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, F_tilde_cols);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), F_tilde_cols, matrix_cols(Ainv_j_inz) - F_tilde_cols);	
			}else{
			    int32_t F_tilde_cols = matrix_cols(step->F_tilde);
			    int32_t Y_cols = matrix_cols(step->Y);

				matrix_free(step->Y);
				step->Y = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), 0, Y_cols);

				matrix_free(step->F_tilde);
				step->F_tilde = matrix_create_sub(Ainv_j_inz, 0, matrix_rows(Ainv_j_inz), Y_cols, matrix_cols(Ainv_j_inz) - Y_cols);
			}
			
			matrix_free(step->R);
			step->R = Ainv_j_j;

		}else{
			kalman_step_index_t inz1 = j/2 - 1;

			kalman_step_index_t physical_column = (kalman_step_index_t) ceil(((double)length)/2) + ptov_n2[inz1];

			matrix_t* L_j_inz = step->F_tilde;

			step_t *physical_step = steps[vtop_n[physical_column]];
			kalman_step_index_t physical_row = physical_step->step;

			matrix_t* Ainv_inz_inz = physical_step->R;

	// 		Ainv_j_inz = -L_j_inz * Ainv_inz_inz;
			matrix_t* Ainv_j_inz = matrix_create(matrix_rows(L_j_inz), matrix_cols(Ainv_inz_inz));
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(-1, L_j_inz, Ainv_inz_inz, 0, Ainv_j_inz);
			
	// 		Dinv = kalman.steps{i}.R;
			matrix_t* Dinv = step->R;

	// 		Ainv_j_j = Dinv * Dinv' - Ainv_j_inz * L_j_inz';
			matrix_t* Dinv_Dinv_t = matrix_create(matrix_rows(Dinv), matrix_rows(Dinv));
			matrix_t * Dinv_t = matrix_create_transpose(Dinv);
			  //printf(">>> %d\n",__LINE__);
			matrix_mutate_gemm(1, Dinv, Dinv_t, 0, Dinv_Dinv_t);
			matrix_free(Dinv_t);
			matrix_t* L_j_inz_t = matrix_create_transpose(L_j_inz);
			matrix_t* tmp = matrix_create(matrix_rows(Ainv_j_inz), matrix_rows(L_j_inz));
			  //printf(">>> %d\n",__LINE__);
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


//void one_layer_converter (void* kalman_v, void* steps_v, int length, int** helper, kalman_step_index_t start, kalman_step_index_t end){
static void one_layer_converter (void* array, kalman_step_index_t length, kalman_step_index_t start, kalman_step_index_t end){
	//kalman_t* kalman = (kalman_t*) kalman_v;
	//step_t* *steps = (step_t**)steps_v;
	
	kalman_step_index_t** helper = (kalman_step_index_t**) array;
	 
	kalman_step_index_t* virtual_to_physical = helper[0];
	kalman_step_index_t* physical_to_virtual = helper[1];
	kalman_step_index_t* variables = helper[2];

	kalman_step_index_t n = variables[0];
	kalman_step_index_t start_ = variables[1];
	kalman_step_index_t jump = variables[2];
	kalman_step_index_t skip = variables[3];
	kalman_step_index_t index = variables[4];

	for (kalman_step_index_t j_ = start; j_ < end; ++j_){
		kalman_step_index_t i = j_ * jump + start_;

		virtual_to_physical[index + j_] = i;
		if (!skip){
			physical_to_virtual[i/2] = (index + j_) - (kalman_step_index_t) ceil(((double)n)/2);
		}
	}
	variables[4] = index;
}

static kalman_step_index_t** virtual_physical(kalman_step_index_t n) {

	kalman_step_index_t* virtual_to_physical = (kalman_step_index_t*) malloc(    n * sizeof(kalman_step_index_t));
	kalman_step_index_t* physical_to_virtual = (kalman_step_index_t*) malloc((n/2) * sizeof(kalman_step_index_t));
	
	kalman_step_index_t* variables = (kalman_step_index_t *)malloc(5 * sizeof(kalman_step_index_t));

	kalman_step_index_t start = 0;
	kalman_step_index_t jump = 2;

	kalman_step_index_t index = 0;

	kalman_step_index_t skip = 1;

	kalman_step_index_t ** result = (kalman_step_index_t **)malloc(3 * sizeof(kalman_step_index_t *));
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
		foreach_in_range(one_layer_converter, result, (kalman_step_index_t) ceil((double)(n - start)/jump), (kalman_step_index_t) ceil((double)(n - start)/jump));

		index += (kalman_step_index_t) ceil((double)(n - start)/jump);

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
static void smooth_recursive(kalman_options_t options, step_t** steps, kalman_step_index_t length) {
	
    if (length == 1) {

		step_t* singleStep = steps[0];

		matrix_t* G = matrix_create_copy(singleStep->G);
		matrix_t* o = matrix_create_copy(singleStep->o);
        //printf(">>> %d\n",__LINE__);

		matrix_t* TAU = matrix_create_mutate_qr(G);
		  //printf(">>> %d\n",__LINE__);
		matrix_mutate_apply_qt(G,TAU,o);

		matrix_mutate_triu(G);

		// ==========================================
		// Cov Change 3
		// ==========================================

		matrix_t* RT = matrix_create_transpose(G);

		matrix_t* RT_R = matrix_create_constant(matrix_cols(G), matrix_cols(G), 0);
		  //printf(">>> %d\n",__LINE__);
		matrix_mutate_gemm(1, RT, G, 0, RT_R);

		singleStep->R = matrix_create_inverse(RT_R);

		singleStep->covariance = matrix_create_copy(singleStep->R);

		matrix_free(RT);
		matrix_free(RT_R);
		
		// ==========================================
		// End Cov Change 3
		// ==========================================
		singleStep->state = matrix_create_trisolve("U",G,o);
		
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
        step_t* step_last = steps[length - 1];
		step_t* step_lmo  = steps[length - 2];

		matrix_t* G_tilde = step_lmo->G_tilde;
		matrix_t* o_i     = step_lmo->c;

		matrix_t* Z       = step_last->Z;
		matrix_t* c_i     = step_last->o;

		matrix_t* QR = matrix_create_vconcat(G_tilde, Z);
        //printf(">>> %d\n",__LINE__);

		matrix_t* TAU = matrix_create_mutate_qr(QR);
		matrix_t* R   = matrix_create_sub(QR,0,matrix_cols(QR), 0, matrix_cols(QR));

		matrix_mutate_triu(R);
		free_and_assign(&(step_lmo->G_tilde), R); //step_lmo->G_tilde = R;

		int32_t ll = MAX(matrix_rows(G_tilde), matrix_cols(G_tilde)); // weird but see Matlab code
		  //printf(">>> %d\n",__LINE__);
		apply_QT_to_block_matrix(QR, TAU, ll, &(step_lmo->c), &(step_last->o));
		matrix_free(QR);
		matrix_free(TAU);
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
	
	kalman_step_index_t** result = virtual_physical(length);

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

#ifdef OBSOLETE
static void smooth(kalman_t* kalman) {

	kalman_step_index_t length = farray_size(kalman->steps);
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
#endif

static void steps_init(void* steps_v, void* steps_array_v, kalman_step_index_t l, kalman_step_index_t start, kalman_step_index_t end) {
  step_t*  steps_array  = (step_t*)   steps_array_v;
  step_t** steps        = (step_t**)  steps_v;
  for (kalman_step_index_t j = start; j < end; j++) {
    steps[j] = &(steps_array[j]);
    step_t* s = steps[j];

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
  }
}

static void steps_weigh(void* equations_v, void* steps_v, kalman_step_index_t l, kalman_step_index_t start, kalman_step_index_t end) {
  kalman_step_equations_t** equations = (kalman_step_equations_t**) equations_v;
  step_t**                  steps     = (step_t**)                  steps_v;
  for (kalman_step_index_t i = start; i < end; i++) {

    steps[i]->step      = equations[i]->step;
    steps[i]->dimension = equations[i]->dimension;

    if (equations[i]->step > 0) {
      matrix_t* H_i = equations[i]->H;
      matrix_t* F_i = equations[i]->F;
      matrix_t* c_i = equations[i]->c;
      matrix_t* K_i = equations[i]->K;
      char K_type   = equations[i]->K_type;

      matrix_t* V_i_H_i = kalman_covariance_matrix_weigh(K_i,K_type,H_i);
      matrix_t* V_i_F_i = kalman_covariance_matrix_weigh(K_i,K_type,F_i);
      matrix_t* V_i_c_i = kalman_covariance_matrix_weigh(K_i,K_type,c_i);

      matrix_mutate_scale(V_i_F_i,-1.0);

      steps[i]->H = V_i_H_i;
      steps[i]->F = V_i_F_i;
      steps[i]->c = V_i_c_i;
    }

    matrix_t* G_i = equations[i]->G;
    matrix_t* o_i = equations[i]->o;
    matrix_t* C_i = equations[i]->C;
    char C_type   = equations[i]->C_type;

    if (o_i != NULL) {
      matrix_t* W_i_G_i = kalman_covariance_matrix_weigh(C_i,C_type,G_i);
      matrix_t* W_i_o_i = kalman_covariance_matrix_weigh(C_i,C_type,o_i);

      steps[i]->G  = W_i_G_i;
      steps[i]->o  = W_i_o_i;
    }
  }
}

static void steps_finalize(void* equations_v, void* steps_v, kalman_step_index_t l, kalman_step_index_t start, kalman_step_index_t end) {
  kalman_step_equations_t** equations = (kalman_step_equations_t**) equations_v;
  step_t**                  steps     = (step_t**)                  steps_v;
  for (kalman_step_index_t i = start; i < end; i++) {

    equations[i]->state      = steps[i]->state;
    equations[i]->covariance = steps[i]->covariance;
    equations[i]->covariance_type = 'C';

    step_t* s = steps[i];

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

    // state and covariance were transfered ...
    //matrix_free(s->state);
    //matrix_free(s->covariance);
  }
}

void kalman_smooth_oddeven(kalman_options_t options, kalman_step_equations_t** equations, kalman_step_index_t l) {
  //step_t** steps = (step_t**) malloc(length * sizeof(step_t*));

  step_t*  steps_array = (step_t*)  malloc( l * sizeof(step_t)  );
  step_t** steps       = (step_t**) malloc( l * sizeof(step_t*) );

  foreach_in_range_two(steps_init,  steps,     steps_array, l, l);
  foreach_in_range_two(steps_weigh, equations, steps,       l, l);

  smooth_recursive(options, steps, l);

  foreach_in_range_two(steps_finalize, equations, steps,       l, l);

  free(steps);
  free(steps_array);
}


#ifdef OBSOLETE
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
#endif

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
