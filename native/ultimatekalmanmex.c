/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier)
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ultimatekalman.h"

#ifndef BUILD_MEX
#error "You must define -DBUILD_MEX to build the mex version of ultimatekalman"
#endif

static double NaN = 0.0 / 0.0;

/*******************************************************************/
/* HANDLES                                                         */
/*******************************************************************/

#define NHANDLES 16
static void* kalman_handles[NHANDLES] = {0 };
static void* cov_handles   [NHANDLES] = { 0 };

static int handleAllocate(void** handles) {
	int i;

	for (i=0; i<NHANDLES; i++) {
		if (handles[i]==NULL) {
			handles[i] = (void*) -1; // mark as being in use with an illegal pointer
			return i;
		}
	}
	return -1;
}

static void handleFree(void** handles, int i) {
	if (i>=0 && i<NHANDLES) handles[i] = NULL;
}

static void handleSet(void** handles, int i, void* p) {
	if (i>=0 && i<NHANDLES) handles[i] = p;
}


static void* handleGet(void** handles, int i) {
	if (i>=0 && i<NHANDLES && handles[i]!=(void*)-1) return handles[i];
	return NULL;
}


/*******************************************************************/
/* MEX HELPER FUNCTIONS                                            */
/*******************************************************************/

static int selector(char* expected, int nrhs, const mxArray *prhs[]) {

  if(nrhs<1) mexErrMsgIdAndTxt("sivantoledo:UltimateKalman:nrhs","mex function requires at least a selector argument.");

  if( !mxIsChar(prhs[0]) ) mexErrMsgIdAndTxt("sivantoledo:UltimateKalman:notChar","selector must be a character array");

  char selector[32];
  //selector = mxGetChars(prhs[0]);

  mxGetString(prhs[0], selector, 32);

  if (expected == NULL) { // a mechanism to display the selector for debugging
  	mexPrintf("selector = <%s>\n",selector);
  	return 0;
  }

  if (!strncmp(selector, expected, strlen(expected))) return 1;

  return 0;
}

static void argCheck(char* selector,
		                 int inMin, int inMax,
										 int outMin, int outMax,
		                 int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[]) {

	int i;
	char msg[81];

  if(nrhs-1<inMin  || nrhs-1>inMax) {
  	sprintf(msg,"Number of inputs to %s must be %d to %d.",selector,inMin,inMax);
  	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",msg);
  }
  if(nlhs<outMin   || nlhs>outMax ) {
  	sprintf(msg,"Number of outputs of %s must be %d to %d.",selector,outMin,outMax);
  	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",msg);
  }

  for (i=1; i<nrhs; i++) if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) { // skip the selector
  	sprintf(msg,"All inputs to %s must be double-precision and real.",selector);
  	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:type",msg);
  }
}

/*******************************************************************/
/* MEX INTERFACES TO ACTUAL FUNCTIONS                              */
/*******************************************************************/

static void mexTemplate(char* selector,int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("create",0,0,1,1,nlhs,plhs,nrhs,prhs);

}


static void mexCreate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("create",0,0,1,1,nlhs,plhs,nrhs,prhs);

  int handle;

  kalman_t* kalman = kalman_create();

  if (kalman==NULL) {
  	handle = -1;
  } else {
  	handle = handleAllocate(kalman_handles);
  	if (handle==-1) kalman_free(kalman);
  	handleSet(kalman_handles,handle,kalman);

  	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  	double* out = mxGetPr(plhs[0]);

  	//printf("create %d %08x\n",handle,kalman);

  	if (handle != -1) out[0] = handle;
  	else              out[0] = NaN;
  }
}

static void mexFree(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("free",1,1,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));

 	void* kalman = handleGet(kalman_handles,handle);
	//printf("free %d %08x\n",handle,kalman);
 	kalman_free( kalman );
 	handleFree(kalman_handles, handle);
}

static void mexEarliest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("earliest",1,1,1,1,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
 	void* p = handleGet(kalman_handles,handle);

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* out = mxGetPr(plhs[0]);

  if (p==NULL) out[0] = NaN;
  else                  out[0] = (double) kalman_earliest(p);
}

static void mexLatest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("latest",1,1,1,1,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
 	void* kalman = handleGet(kalman_handles,handle);

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* out = mxGetPr(plhs[0]);

  if (kalman==NULL) out[0] = NaN;
  else                  out[0] = (double) kalman_latest(kalman);
}

static void mexEvolve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("evolve",7,7,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	int32_t n_i = (int32_t) mxGetScalar(prhs[2]);

	//printf("evolve %d %d %08x\n",n_i,handle,kalman);

	matrix_t* H_i = matrix_create_from_mxarray(prhs[3]);
	matrix_t* F_i = matrix_create_from_mxarray(prhs[4]);
	matrix_t* c_i = matrix_create_from_mxarray(prhs[5]);
	matrix_t* K_i = matrix_create_from_mxarray(prhs[6]);

	char type = (char) round(mxGetScalar(prhs[7]));

	kalman_evolve(kalman, n_i, H_i, F_i, c_i, K_i, type);

	//return;

	matrix_free(K_i);
	matrix_free(H_i);
	matrix_free(F_i);
	matrix_free(c_i);
}

static void mexObserve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("evolve",5,5,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	matrix_t* G_i = matrix_create_from_mxarray(prhs[2]);
	matrix_t* o_i = matrix_create_from_mxarray(prhs[3]);
	matrix_t* C_i = matrix_create_from_mxarray(prhs[4]);

	char type = (char) round(mxGetScalar(prhs[5]));

	kalman_observe(kalman, G_i, o_i, C_i, type);

	matrix_free(C_i);
	matrix_free(o_i);
	matrix_free(G_i);
}

static void mexSmooth(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("smooth",1,1,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));

 	void* kalman = handleGet(kalman_handles,handle);
	//printf("smooth %d %08x\n",handle,kalman);
 	//kalman_free( kalman );
 	//handleFree(kalman_handles, handle);
}

static void mexEstimate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("estimate",2,2,1,1,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	int64_t s_i = (int64_t) mxGetScalar(prhs[2]);

	//printf("estimate %d %d %08x\n",s_i,handle,kalman);

	matrix_t* e = kalman_estimate(kalman,s_i);

	if (e == NULL) {
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
	} else {
		plhs[0] =  matrix_copy_to_mxarray(e);
		matrix_free(e);
	}
}

static void mexCovariance(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("covariance",2,2,1,1,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	int64_t s_i = (int64_t) mxGetScalar(prhs[2]);

	//printf("covariance %d %d %08x\n",s_i,handle,kalman);

	matrix_t* W = kalman_covariance(kalman,s_i);

	if (W == NULL) {
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
	} else {
		plhs[0] =  matrix_copy_to_mxarray(W);
		matrix_free(W);
	}
}

static void mexForget(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("forget",2,2,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	int64_t s_i = (int64_t) mxGetScalar(prhs[2]);

	printf("forget %d %d %08x\n",s_i,handle,kalman);

	kalman_forget(kalman,s_i);
}

static void mexRollback(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	argCheck("rollback",2,2,0,0,nlhs,plhs,nrhs,prhs);

	int handle = (int) floor(mxGetScalar(prhs[1]));
	void* kalman = handleGet(kalman_handles,handle);

	int64_t s_i = (int64_t) mxGetScalar(prhs[2]);

	printf("rollback %d %d %08x\n",s_i,handle,kalman);

	kalman_rollback(kalman,s_i);
}



/*******************************************************************/
/* MEX FUNCTION                                                    */
/*******************************************************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    //mexPrintf("Build time = %s %s\n",__DATE__,__TIME__);


    selector(NULL,nrhs,prhs);


    if      (selector("create"     ,nrhs,prhs)) { mexCreate    (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("free"       ,nrhs,prhs)) { mexFree      (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("earliest"   ,nrhs,prhs)) { mexEarliest  (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("latest"     ,nrhs,prhs)) { mexLatest    (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("evolve"     ,nrhs,prhs)) { mexEvolve    (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("observe"    ,nrhs,prhs)) { mexObserve   (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("smooth"     ,nrhs,prhs)) { mexSmooth    (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("estimate"   ,nrhs,prhs)) { mexEstimate  (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("covariance" ,nrhs,prhs)) { mexCovariance(nlhs,plhs,nrhs,prhs); return; }
    else if (selector("forget"     ,nrhs,prhs)) { mexForget    (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("rollback"   ,nrhs,prhs)) { mexRollback  (nlhs,plhs,nrhs,prhs); return; }
    else mexErrMsgIdAndTxt("sivantoledo:UltimateKalman:invalidSelector","invalid selector");
}
