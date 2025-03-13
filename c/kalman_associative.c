
/*
 * kalman_associative.c
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

#define KALMAN_MATRIX_SHORT_TYPE
#include "kalman.h"
#include "parallel.h"

#ifdef PARALLEL

#ifdef MACOS
void* local_aligned_alloc(size_t alignment, size_t size) {
	void* p;
	int result = posix_memalign(&p, alignment, size);
	if (result != 0) return NULL;
	return p;
} 
#define malloc(x) local_aligned_alloc(64,(x))
#else
#define malloc(x) aligned_alloc(64,(x))
#endif

//#define parallel_for_c parallel_for_c_associative
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
	
	char   C_type;
	char   K_type;

	kalman_matrix_t* H;
	kalman_matrix_t* F;
	kalman_matrix_t* K;
	kalman_matrix_t* c;

	kalman_matrix_t* G;
	kalman_matrix_t* o;
	kalman_matrix_t* C;

	kalman_matrix_t* Z;

	kalman_matrix_t* A;
	kalman_matrix_t* b;

	kalman_matrix_t* e;
	kalman_matrix_t* J;

	kalman_matrix_t* E;
	kalman_matrix_t* g;
	kalman_matrix_t* L;

	kalman_matrix_t* state;
	kalman_matrix_t* covariance;
} step_t;

static void* step_create() {
	step_t* s = malloc(sizeof(step_t));
	s->step      = -1;
	s->dimension = -1;
	
	s->C_type    = 0;
	s->K_type    = 0;


	s->H = NULL;
	s->F = NULL;
	s->K = NULL;
	s->c = NULL;

	s->G = NULL;
	s->o = NULL;
	s->C = NULL;

	s->Z = NULL;

	s->A = NULL;
	s->b = NULL;

	s->e = NULL;
	s->J = NULL;

	s->E = NULL;
	s->g = NULL;
	s->L = NULL;
	
	s->state = NULL;
	s->covariance = NULL;

	assert( s != NULL );
	return s;
}

static void step_free(void* v) {
	step_t* s = (step_t*) v;
	
	matrix_free(s->H);
	matrix_free(s->F);
	matrix_free(s->K);
	matrix_free(s->c);

	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->Z);

	matrix_free(s->A);
	matrix_free(s->b);

	matrix_free(s->e);
	matrix_free(s->J);

	matrix_free(s->E);
	matrix_free(s->g);
	matrix_free(s->L);
	
	matrix_free(s->state);
	matrix_free(s->covariance);

	free(s);
}

static void step_rollback(void* v) {
	step_t* s = (step_t*) v;

	matrix_free(s->G);
	matrix_free(s->o);
	matrix_free(s->C);

	matrix_free(s->Z);

	matrix_free(s->A);
	matrix_free(s->b);

	matrix_free(s->e);
	matrix_free(s->J);

	matrix_free(s->E);
	matrix_free(s->g);
	matrix_free(s->L);
	
	matrix_free(s->state);
	matrix_free(s->covariance);
}

static int64_t step_get_index(void* v) {
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


	// matrix_t* V_i_H_i = cov_weigh(K_i,K_type,H_i);
	// matrix_t* V_i_F_i = cov_weigh(K_i,K_type,F_i);
	// matrix_t* V_i_c_i = cov_weigh(K_i,K_type,c_i);

	// matrix_mutate_scale(V_i_F_i,-1.0);

	kalman_current->H	= matrix_create_copy(H_i);
	kalman_current->F 	= matrix_create_copy(F_i);
	kalman_current->c  = matrix_create_copy(c_i);
	kalman_current->K  = matrix_create_copy(K_i);
	kalman_current->K_type = K_type;
	
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
	printf("observe %d\n",(int) kalman->current->step);
#endif

	// matrix_t* W_i_G_i = NULL;
	// matrix_t* W_i_o_i = NULL;

	if (o_i != NULL) {
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

	kalman_current->G  = matrix_create_copy(G_i);
	kalman_current->o  = matrix_create_copy(o_i);
	kalman_current->C  = matrix_create_copy(C_i);
	kalman_current->C_type = C_type;

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
/* CONCURRENT SET OF POINTERS                                                 */
/******************************************************************************/
/*
 * We define a simple concurrent set data structure, to keep track of steps
 * that are created during the parallel prefix sum operation (the elements
 * are structures, not value types) so that we can release them at the end of
 * the operation.
 */
 

/*
 * FNV-1a hash of an address, to generate a random integer
 */

uint32_t hash_uint32(uint32_t value) {
    const uint32_t FNV_PRIME = 16777619u;
    const uint32_t FNV_OFFSET = 3141592653u; // Different offset for distinctness

    uint32_t hash = FNV_OFFSET;
    for (size_t i = 0; i < sizeof(uint32_t); i++) {
        hash ^= (value & 0xFF);
        hash *= FNV_PRIME;
        value >>= 8;
    }
    return hash;
}

typedef struct concurrent_set_st {
  int            size;
  void**         pointers;          
  spin_mutex_t** locks; // pointers to locks, to allow locks of any type
  void   (*foreach)(void*);
} concurrent_set_t;

//static void parallelInit(void* la_v, void* *helper, size_t length, size_t start, size_t end){
static void concurrent_set_parallel_init(void* set_v, int length, size_t start, size_t end){
    concurrent_set_t* set = (concurrent_set_t*) set_v;
	for (int i = start; i < end; i++) {
        (set->pointers)[i] = NULL;
        (set->locks)[i] = spin_mutex_create();
    }
}

static void concurrent_set_parallel_destroy(void* set_v, int length, size_t start, size_t end){
    concurrent_set_t* set = (concurrent_set_t*) set_v;
	for (int i = start; i < end; i++) {
        spin_mutex_destroy((set->locks)[i]);
    }
}

static void concurrent_set_parallel_foreach(void* set_v, int length, size_t start, size_t end){
    concurrent_set_t* set = (concurrent_set_t*) set_v;
	for (int i = start; i < end; i++) {
	  if ( (set->pointers)[i] != NULL ) {
	    (*(set->foreach))((set->pointers)[i]);
	  }
    }
}

static concurrent_set_t* concurrent_set_create(int capacity, void (*foreach)(void*)) {
    concurrent_set_t* set = (concurrent_set_t*) malloc(sizeof(concurrent_set_t));
    set->size = capacity * 10; // expansion to reduce contention
    set->foreach = foreach;
    set->pointers = (void**)         malloc( (set->size) * sizeof(void*));
    set->locks    = (spin_mutex_t**) malloc( (set->size) * sizeof(spin_mutex_t*));

    //parallel_for_c(la, NULL, 0, k, BLOCKSIZE, parallelInit);
    foreach_in_range(concurrent_set_parallel_init, set, set->size, set->size);

    return set;
}

static void concurrent_set_free(concurrent_set_t* set) {
  foreach_in_range(concurrent_set_parallel_destroy, set, set->size, set->size);
  //parallel_for_c(la, NULL, 0, la->rows, BLOCKSIZE, parallelDestroy);
  free(set->locks);
  free(set->pointers);
  free(set);
}

static void concurrent_set_insert(concurrent_set_t* set, void* element) {
	uint32_t inserted = 0;
	uint32_t h = (uintptr_t) element;
	uint32_t i;
	do {
		h = hash_uint32(h);
		i = h % (set->size);
		
	  	spin_mutex_lock((set->locks)[i]);
	  	if ((set->pointers)[i] == NULL) {
			(set->pointers)[i] = element;
			inserted = 1;
		}
	  	spin_mutex_unlock((set->locks)[i]);
	} while (!inserted);
}

static void concurrent_set_foreach(concurrent_set_t* set) {
  foreach_in_range(concurrent_set_parallel_foreach, set, set->size, set->size);
}

/******************************************************************************/
/* Associative Smoother                                                       */
/******************************************************************************/

static void build_filtering_element(kalman_t* kalman, int i) {
	/*
		we denote by Z the matrix denoted by C in the article,
		because we already use C for the covariance of the
		observations.
	*/

	step_t* step_i = farray_get(kalman->steps, i);
	int n_i = step_i->dimension;
	int step = step_i->step;


	if (step == 0) {
		return;
	}

	if (step == 1) {
		step_t* step_1 = farray_get(kalman->steps, 0);
		matrix_t* G_i = step_1->G;
		matrix_t* o_i = step_1->o;
		matrix_t* C_i = step_1->C;

		char C_type = step_i->C_type;
						
		matrix_t* W_i_G_i = cov_weigh(C_i, C_type, G_i);
		matrix_t* W_i_o_i = cov_weigh(C_i, C_type, o_i);

		
		matrix_t* R = matrix_create_copy(W_i_G_i);
		matrix_t* Q = matrix_create_mutate_qr(R);
		
		matrix_mutate_apply_qt(R,Q,W_i_o_i);
		
		matrix_mutate_triu(R);
		matrix_t* m0 = matrix_create_trisolve(R,W_i_o_i);

		matrix_t* RT = matrix_create_transpose(R);
		matrix_t* RTR = matrix_create_multiply(RT,R);
		matrix_t* P0 = matrix_create_inverse(RTR);


		step_1->state = m0;
		step_1->covariance = P0;
		
		matrix_free(Q);
		matrix_free(R);

		matrix_free(W_i_G_i);
		matrix_free(W_i_o_i);
		
		matrix_free(RT);
		matrix_free(RTR);
	}

	matrix_t* F_i = step_i->F;
	matrix_t* c_i = step_i->c;
	matrix_t* K_i = explicit(step_i->K,step_i->K_type);

	if (step == 1) {
		step_t* step_1 = farray_get(kalman->steps, 0);
		matrix_t* F_iT = matrix_create_transpose(F_i);
		matrix_t* P0 = step_1->covariance;
		matrix_t* FP0 = matrix_create_multiply(F_i,P0);
		matrix_t* FPFT = matrix_create_multiply(FP0,F_iT);
		matrix_t* new_K_i = matrix_create_add(K_i,FPFT);
		matrix_free(K_i);
		K_i = new_K_i;
		matrix_free(F_iT);
		matrix_free(FP0);
		matrix_free(FPFT);
	}

	if (step_i->o == NULL) {
		step_i->Z = K_i;
		if (step == 1) {
			step_i->A = matrix_create_constant(n_i,n_i,0.0);
			step_t* step_1 = farray_get(kalman->steps, 0);
			matrix_t* m0 = step_1->state;
			matrix_t* b = matrix_create_add(m0,c_i);
			step_i->b = b;
		}else{
			step_i->A = matrix_create_copy(F_i);	
			step_i->b = matrix_create_copy(c_i);
		}
		step_i->e = NULL;
		step_i->J = NULL;
	}else{ // there are observations
		matrix_t* G_i = step_i->G;
		matrix_t* o_i = step_i->o;
		matrix_t* C_i = explicit(step_i->C,step_i->C_type);

		matrix_t* G_iT = matrix_create_transpose(G_i);
		matrix_t* KGT = matrix_create_multiply(K_i,G_iT);
		matrix_t* GKGT = matrix_create_multiply(G_i,KGT);
		matrix_t* S = matrix_create_add(GKGT,C_i);

		matrix_free(G_iT);
		matrix_free(KGT);	
		matrix_free(GKGT);
		matrix_free(C_i);

		matrix_t* ST = matrix_create_transpose(S);

		matrix_t* G_i_trans_inv_S_T = matrix_create_mldivide(ST,G_i); 
		matrix_t* G_i_trans_inv_S = matrix_create_transpose(G_i_trans_inv_S_T);
		
		matrix_free(ST);
		matrix_free(G_i_trans_inv_S_T);
		
		matrix_t* K = matrix_create_multiply(K_i, G_i_trans_inv_S);

		if (step == 1) {
			step_i->A = matrix_create_constant(n_i,n_i,0.0);
			step_t* step_1 = farray_get(kalman->steps, 0);
			matrix_t* m0 = step_1->state;
			matrix_t* F_im = matrix_create_multiply(F_i,m0);
			matrix_t* m1 = matrix_create_add(F_im,c_i);
			matrix_t* G_im = matrix_create_multiply(G_i,m1);
			matrix_t* o_G_im = matrix_create_subtract(o_i,G_im);
			matrix_t* K_o_G_im = matrix_create_multiply(K,o_G_im);
			matrix_t* b = matrix_create_add(m1,K_o_G_im);
			step_i->b = b;

			matrix_t *KS = matrix_create_multiply(K,S);
			matrix_t* KT = matrix_create_transpose(K);
			matrix_t* KSKT = matrix_create_multiply(KS,KT);
			matrix_t* Z = matrix_create_subtract(K_i,KSKT);
			step_i->Z = Z;

			matrix_free(F_im);
			matrix_free(m1);
			matrix_free(G_im);
			matrix_free(o_G_im);
			matrix_free(K_o_G_im);

			matrix_free(KS);
			matrix_free(KT);
			matrix_free(KSKT);
		}else{
			matrix_t* GF = matrix_create_multiply(G_i,F_i);
			matrix_t* KGF = matrix_create_multiply(K,GF);
			matrix_t* A = matrix_create_subtract(F_i,KGF);
			step_i->A = A;

			matrix_free(GF);
			matrix_free(KGF);

			matrix_t* G_ic = matrix_create_multiply(G_i,c_i);
			matrix_t* o_G_ic = matrix_create_subtract(o_i,G_ic);
			matrix_t* K_o_G_ic = matrix_create_multiply(K,o_G_ic);
			matrix_t* b = matrix_create_add(c_i,K_o_G_ic);
			step_i->b = b;

			matrix_free(G_ic);
			matrix_free(o_G_ic);
			matrix_free(K_o_G_ic);

			matrix_t* KG = matrix_create_multiply(K,G_i);
			matrix_t* KGK_i = matrix_create_multiply(KG,K_i);
			matrix_t* Z = matrix_create_subtract(K_i,KGK_i);
			step_i->Z = Z;

			matrix_free(KG);
			matrix_free(KGK_i);
		}
		matrix_free(K_i);

		matrix_t* G_ic = matrix_create_multiply(G_i,c_i);
		matrix_t* o_G_ic = matrix_create_subtract(o_i,G_ic);
		matrix_t* FT = matrix_create_transpose(F_i);
		matrix_t* FTG = matrix_create_multiply(FT,G_i_trans_inv_S);
		matrix_t* e = matrix_create_multiply(FTG,o_G_ic);
		step_i->e = e;

		matrix_free(G_ic);
		matrix_free(o_G_ic);
		matrix_free(FT);
		
		matrix_t* GF = matrix_create_multiply(G_i,F_i);
		matrix_t* J = matrix_create_multiply(FTG,GF);
		step_i->J = J;

		matrix_free(FTG);
		matrix_free(GF);

		matrix_free(G_i_trans_inv_S);
		matrix_free(K);
		matrix_free(S);
	}
}

static
void build_smoothing_element(kalman_t* kalman, int i) {

	if (i == farray_size(kalman->steps) - 1) {
		step_t* step_i = farray_get(kalman->steps,i);
		int ni = step_i->dimension;
		step_i->E = matrix_create_constant(ni,ni,0.0);
		step_i->g = matrix_create_copy(step_i->state);
		step_i->L = matrix_create_copy(step_i->covariance);
	}else{
		step_t* step_i = farray_get(kalman->steps,i);
		matrix_t* x = step_i->state;
		matrix_t* P = explicit(step_i->covariance,'C');
		step_t* step_ip1 = farray_get(kalman->steps,i+1);
		matrix_t* F = step_ip1->F;
		matrix_t* Q = explicit(step_ip1->K,step_ip1->K_type);
		matrix_t* c = step_ip1->c;

		matrix_t* FT = matrix_create_transpose(F);
		matrix_t* PFT = matrix_create_multiply(P,FT);
		matrix_t* FPFT = matrix_create_multiply(F,PFT);
		matrix_t* FPFT_Q = matrix_create_add(FPFT,Q);

		matrix_t* PFT_T = matrix_create_transpose(PFT); 
		matrix_t* FPFT_Q_T = matrix_create_transpose(FPFT_Q);
		matrix_t* E_T = matrix_create_mldivide(FPFT_Q_T,PFT_T);
		matrix_t* E = matrix_create_transpose(E_T);

		step_i->E = E;

		matrix_free(FT);
		matrix_free(PFT);
		matrix_free(FPFT);
		matrix_free(FPFT_Q);
		matrix_free(PFT_T);
		matrix_free(FPFT_Q_T);
		matrix_free(E_T);


		matrix_t* Fx = matrix_create_multiply(F,x);
		matrix_t* Fx_c = matrix_create_add(Fx,c);
		matrix_t* E_Fx_c = matrix_create_multiply(E,Fx_c);
		matrix_t* g = matrix_create_subtract(x,E_Fx_c);
		step_i->g = g;

		matrix_free(Fx);
		matrix_free(Fx_c);
		matrix_free(E_Fx_c);
		matrix_t* EF = matrix_create_multiply(E,F);
		matrix_t* EFP = matrix_create_multiply(EF,P);
		matrix_t* L = matrix_create_subtract(P,EFP);
		step_i->L = L;

		matrix_free(EF);
		matrix_free(EFP);
		
		matrix_free(P);
		matrix_free(Q);
	}
}

static void build_filtering_elements(void* kalman_v, int l, size_t start, size_t end) {
    kalman_t* kalman = (kalman_t*) kalman_v;
    
	for (int j=start; j < end; j++) {
		build_filtering_element(kalman,j);
	}
}

static void build_smoothing_elements(void* kalman_v, int l, size_t start, size_t end) {
    kalman_t* kalman = (kalman_t*) kalman_v;
	
	for (int j=start; j < end; j++) {
		build_smoothing_element(kalman,j);
	}
}

//step_t* filteringAssociativeOperation(step_t* si, step_t* sj, concurrent_set_t* created_steps, int row, int is_final_scan) {
static void* filteringAssociativeOperation(void* si_v, void* sj_v, void* created_steps_v, int is_final_scan) {
  step_t* si = (step_t*) si_v;
  step_t* sj = (step_t*) sj_v;
  
  if (si == NULL){
    return sj;
  }
  
  if (sj == NULL){
    return si;
  }


  step_t* sij = step_create();

#ifdef PARALLEL
  concurrent_set_t* created_steps = (concurrent_set_t*) created_steps_v;
  if (!is_final_scan){
    concurrent_set_insert(created_steps, sij);
  }
#endif
		int ni = matrix_rows(si->b);

		matrix_t* eye_ni = matrix_create_identity(ni,ni);
		matrix_t* siZ_sjJ = matrix_create_multiply(si->Z,sj->J);
		matrix_t* eye_ni_plus_siZ_sjJ = matrix_create_add(eye_ni,siZ_sjJ);

		matrix_t* AT = matrix_create_transpose(sj->A);
		matrix_t* other_T = matrix_create_transpose(eye_ni_plus_siZ_sjJ);
		matrix_t* XT = matrix_create_mldivide(other_T,AT);
		matrix_t* X = matrix_create_transpose(XT);
		
		
		matrix_free(siZ_sjJ);
		matrix_free(eye_ni_plus_siZ_sjJ);
		matrix_free(AT);
		matrix_free(other_T);
		matrix_free(XT);


		matrix_t* sjJ_siZ = matrix_create_multiply(sj->J,si->Z);
		matrix_t* eye_ni_plus_sjJ_siZ = matrix_create_add(eye_ni,sjJ_siZ);
		matrix_t* A_iT = matrix_create_transpose(si->A);

		
		other_T = matrix_create_transpose(eye_ni_plus_sjJ_siZ);
		matrix_t* YT = matrix_create_mldivide(other_T,si->A);
		matrix_t* Y = matrix_create_transpose(YT);

		matrix_free(sjJ_siZ);
		matrix_free(eye_ni_plus_sjJ_siZ);
		matrix_free(eye_ni);
		matrix_free(A_iT);
		matrix_free(other_T);
		matrix_free(YT);

		matrix_t* XA = matrix_create_multiply(X,si->A);
		sij->A = XA;

		matrix_t* siZ_sj_e = matrix_create_multiply(si->Z,sj->e);
		matrix_t* siZ_sj_e_plus_si_b = matrix_create_add(siZ_sj_e,si->b);
		matrix_t* X_siZ_sj_e_plus_si_b = matrix_create_multiply(X,siZ_sj_e_plus_si_b);
		matrix_t* b = matrix_create_add(X_siZ_sj_e_plus_si_b,sj->b);
		sij->b = b;

		matrix_free(siZ_sj_e);
		matrix_free(siZ_sj_e_plus_si_b);
		matrix_free(X_siZ_sj_e_plus_si_b);

		matrix_t* X_siZ = matrix_create_multiply(X,si->Z);
		matrix_t* A_jT = matrix_create_transpose(sj->A);
		matrix_t* X_siZ_AT = matrix_create_multiply(X_siZ,A_jT);
		matrix_t* Z = matrix_create_add(X_siZ_AT,sj->Z);
		sij->Z = Z;

		matrix_free(X_siZ);
		matrix_free(A_jT);
		matrix_free(X_siZ_AT);

		matrix_t* sjJ_si_b = matrix_create_multiply(sj->J,si->b);
		matrix_t* sj_e_minus_sjJ_si_b = matrix_create_subtract(sj->e,sjJ_si_b);
		matrix_t* Y_sj_e_minus_sjJ_si_b = matrix_create_multiply(Y,sj_e_minus_sjJ_si_b);
		matrix_t* e = matrix_create_add(Y_sj_e_minus_sjJ_si_b,si->e);
		sij->e = e;

		matrix_free(sjJ_si_b);
		matrix_free(sj_e_minus_sjJ_si_b);
		matrix_free(Y_sj_e_minus_sjJ_si_b);

		matrix_t* sjJ_siA = matrix_create_multiply(sj->J,si->A);
		matrix_t* Y_sjJ_siA = matrix_create_multiply(Y,sjJ_siA);
		matrix_t* J = matrix_create_add(Y_sjJ_siA,si->J);
		sij->J = J;

		matrix_free(sjJ_siA);
		matrix_free(Y_sjJ_siA);
		matrix_free(X);
		matrix_free(Y);

		return sij;

}


//step_t* smoothingAssociativeOperation(step_t* si, step_t* sj, concurrent_set_t* created_steps, int row, int is_final_scan) {
static void* smoothingAssociativeOperation(void* si_v, void* sj_v, void* created_steps_v, int is_final_scan) {
  step_t* si = (step_t*) si_v;
  step_t* sj = (step_t*) sj_v;

  if (si == NULL){
    return sj;
  }
  
  if (sj == NULL){
    return si;
  }
  
  step_t* sij = step_create();
#ifdef PARALLEL
  concurrent_set_t* created_steps = (concurrent_set_t*) created_steps_v;

  if (!is_final_scan){
    concurrent_set_insert(created_steps, sij);
  }
#endif
		
		matrix_t* E = matrix_create_multiply(sj->E,si->E);
		sij->E = E;

		matrix_t* Eg = matrix_create_multiply(sj->E,si->g);
		matrix_t* g = matrix_create_add(Eg,sj->g);
		sij->g = g;

		matrix_free(Eg);

		matrix_t* ET = matrix_create_transpose(sj->E);
		matrix_t* EL = matrix_create_multiply(sj->E,si->L);
		matrix_t* ELT = matrix_create_multiply(EL,ET);
		matrix_t* L = matrix_create_add(ELT,sj->L);
		sij->L = L;

		matrix_free(ET);
		matrix_free(EL);
		matrix_free(ELT);

		return sij;
}


static void filtered_to_state(void* kalman_v, void* filtered_v, size_t l, size_t start, size_t end){
  kalman_t* kalman   = (kalman_t*) kalman_v;
  step_t**  filtered = (step_t**)  filtered_v;

	int i = start;
	for (int j = start + 1; j < end + 1; j++) {
		step_t* step_j = farray_get(kalman->steps,j);
		step_t* filtered_i = filtered[i];
		step_j->state = matrix_create_copy(filtered_i->b);
		step_j->covariance = matrix_create_copy(filtered_i->Z);		
		if (i != 0) {
			step_free(filtered_i);
		}
		i++;
	}
}

static void smoothed_to_state(void* kalman_v, void* smoothed_v, size_t l, size_t start, size_t end){
  kalman_t* kalman   = (kalman_t*) kalman_v;
  step_t**  smoothed = (step_t**)  smoothed_v;
  
	int i = 0;
	for (int j = start; j < end; j++) {
		i = l - 2 - j + 1;
		step_t* step_j = farray_get(kalman->steps,j);
		step_t* smoothed_i = smoothed[i];
		matrix_free(step_j->state);
		step_j->state = matrix_create_copy(smoothed_i->g);
		matrix_free(step_j->covariance);
		step_j->covariance = matrix_create_copy(smoothed_i->L);
		step_free(smoothed_i);
	}
}

//step_t** cummulativeSumsSequential(kalman_t* kalman, step_t* (*f)(step_t*, step_t*), int s, int e, int stride) {
static step_t** cummulativeSumsSequential(kalman_t* kalman, void* (*f)(void*, void*, void*, int, int), int s, int e, int stride) {
	step_t** sums = (step_t**)malloc((abs(e-s) + 1)*sizeof(step_t*));
	int i = 0;
	step_t** a = (step_t**)kalman->steps->elements;
	step_t* sum = a[s];
	sums[i] = sum;
	i++;

	if (s > e) {
		for (int j=s+stride; j>=e; j+=stride) {
		        sum = (step_t*) f(sum,a[j], NULL, -1, -1);
			sums[i] = sum;
			i++;
		}
	} else {
		for (int j=s+stride; j<=e; j+=stride) {
			sum = (step_t*) f(sum,a[j], NULL, -1, -1);
			sums[i] = sum;
			i++;
		}
	}
	return sums;
}

static void prefix_sums_sequential(void* (*f)(void*, void*, void*, int, int), void** a, void** sums, int s, int e, int stride) {
	//step_t** sums = (step_t**)malloc((abs(e-s) + 1)*sizeof(step_t*));
	int i = 0;
	//step_t** a = (step_t**)kalman->steps->elements;
	//step_t* sum = a[s];
	void* sum = a[s];
	sums[i] = sum;
	i++;

	if (s > e) {
		for (int j=s+stride; j>=e; j+=stride) {
		    sum = f(sum,a[j], NULL, -1, -1);
			sums[i] = sum;
			i++;
		}
	} else {
		for (int j=s+stride; j<=e; j+=stride) {
			sum = f(sum,a[j], NULL, -1, -1);
			sums[i] = sum;
			i++;
		}
	}
	return sums;
}


static void smooth(kalman_t* kalman) {
	int l = farray_size(kalman->steps);

	
//#ifdef PARALLEL
//	parallel_for_c(kalman, NULL, l, l, BLOCKSIZE, buildFilteringElements);
//#else
//	buildFilteringElements(kalman, NULL, l, 0, l);
//#endif
	foreach_in_range(build_filtering_elements, kalman, l, l);



	step_t** filtered = (step_t**) malloc((l-1) * sizeof(step_t*));
#ifdef PARALLEL
	concurrent_set_t* filtered_created_steps = concurrent_set_create(l, step_free);
	parallel_scan_c(kalman->steps->elements, (void**) filtered, filtered_created_steps, filteringAssociativeOperation , l - 1, 1);
	concurrent_set_foreach(filtered_created_steps);
	concurrent_set_free(filtered_created_steps);
#else
	//step_t** filtered = cummulativeSumsSequential(kalman, filteringAssociativeOperation, 1, l - 1, 1);
	prefix_sums_sequential(filteringAssociativeOperation, kalman->steps->elements, filtered, 1, l-1, 1);
#endif

//#ifdef PARALLEL
//	parallel_for_c(kalman, (void**) filtered, l, l - 1, BLOCKSIZE, filtered_to_state);
//#else
//	filtered_to_state(kalman, (void**) filtered, l, 0, l - 1);
//#endif
	foreach_in_range_two(filtered_to_state, kalman, filtered, l, l-1);

	free(filtered);

//#ifdef PARALLEL
//	parallel_for_c(kalman, NULL, l, l, BLOCKSIZE, buildSmoothingElements);
//#else
//	buildSmoothingElements(kalman, NULL, l, 0, l);
//#endif
	foreach_in_range(build_smoothing_elements, kalman, l, l);


	step_t** smoothed = (step_t**) malloc(l * sizeof(step_t*));
#ifdef PARALLEL
	concurrent_set_t* smoothed_created_steps = concurrent_set_create(l, step_free);
	parallel_scan_c(kalman->steps->elements, (void**) smoothed, smoothed_created_steps, smoothingAssociativeOperation , l, -1);
	concurrent_set_foreach(smoothed_created_steps);
	concurrent_set_free(smoothed_created_steps);
#else
	//step_t** smoothed = cummulativeSumsSequential(kalman, smoothingAssociativeOperation, l - 1, 0, -1);
	prefix_sums_sequential(smoothingAssociativeOperation, kalman->steps->elements, smoothed, l-1, 0, -1);
#endif

//#ifdef PARALLEL
//	parallel_for_c(kalman, (void**) smoothed, l, l - 1, BLOCKSIZE, smoothed_to_state);
//#else
//	smoothed_to_state(kalman, (void**) smoothed, l, 0, l - 1);
//#endif
	foreach_in_range_two(smoothed_to_state, kalman, smoothed, l, l-1);

	free(smoothed);
}

void kalman_create_associative(kalman_t* kalman) {
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
