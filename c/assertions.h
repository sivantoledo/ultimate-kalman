#ifndef ASSERTIONS_H
#define ASSERTIONS_H

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

#endif
