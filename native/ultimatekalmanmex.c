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

static double NaN = 0.0 / 0.0;

/* The computational routine */
void arrayProduct(double x, double *y, double *z, mwSize n)
{
    mwSize i;
    /* multiply each element y by x */
    for (i=0; i<n; i++) {
        z[i] = x * y[i];
    }
}

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

  for (i=1; i<nrhs; i++) if( !mxIsDouble(prhs[i]) ) { // skip the selector
  	sprintf(msg,"All inputs to %s must be real.",selector);
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

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double* out = mxGetPr(plhs[0]);

  int handle = handleAllocate(kalman_handles);

  handleSet(kalman_handles,handle,kalman_create());

  if (handle != -1) out[0] = handle;
  else              out[0] = NaN;
}

static void mexRelease(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //mexPrintf("-- release\n");
	argCheck("release",1,1,0,0,nlhs,plhs,nrhs,prhs);

   //int handle = (int) floor(phandle[0]);

   int handle = (int) floor(mxGetScalar(prhs[1]));

   handleFree(kalman_handles, handle);

   //mexPrintf("-- release handle=%d\n",handle);
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

    mexPrintf("Build time = %s %s\n",__DATE__,__TIME__);

    selector(NULL,nrhs,prhs);


    if      (selector("create"   ,nrhs,prhs)) { mexCreate  (nlhs,plhs,nrhs,prhs); return; }
    else if (selector("release"  ,nrhs,prhs)) { mexRelease (nlhs,plhs,nrhs,prhs); return; }
    /*
    else if (!strcmp(selector,"earliest")) {
      mexPrintf("-- earliest\n");
    }

    else if (!strcmp(selector,"latest  ")) {
      mexPrintf("-- latest\n");
    }

    else if (!strcmp(selector,"evolve  ")) {
      mexPrintf("-- evolve\n");
    }

    else if (!strcmp(selector,"observe ")) {
      mexPrintf("-- observe\n");
    }

    else if (!strcmp(selector,"estimate")) {
      mexPrintf("-- estimate\n");
    }

    else if (!strcmp(selector,"smooth  ")) {
      mexPrintf("-- smooth\n");
    }

    else if (!strcmp(selector,"forget  ")) {
      mexPrintf("-- forget\n");
    }

    else if (!strcmp(selector,"rollback")) {
      mexPrintf("-- rollback\n");
    }
    */
    else mexErrMsgIdAndTxt("sivantoledo:UltimateKalman:invalidSelector","invalid selector");

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }

    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }

    /* get the value of the scalar input  */
    multiplier = mxGetScalar(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[1]);
    #else
    inMatrix = mxGetPr(prhs[1]);
    #endif

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
    #else
    outMatrix = mxGetPr(plhs[0]);
    #endif

    /* call the computational routine */
    arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);
}
