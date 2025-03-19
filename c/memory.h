#ifndef MEMORY_H
#define MEMORY_H

#if defined(BUILD_MEX) && defined(BUILD_MATLAB)
#include <stdlib.h>

#else // not MATLAB MEX

#include <stdlib.h>

// Helper macro to pad size to the next multiple of 64
#define PAD_TO_64(x) (((x) + 63) & ~63)

#if defined(_WIN32) || defined(_WIN64) // Windows
    #include <malloc.h>
    #define malloc(x) _aligned_malloc(PAD_TO_64(x), 64)
#elif defined(__APPLE__) // macOS
    #define malloc(x) ({ \
        void *ptr; \
        posix_memalign(&ptr, 64, PAD_TO_64(x)); \
        ptr; \
    })
#else // Linux (and other POSIX systems with C11)
    #define malloc(x) aligned_alloc(64, PAD_TO_64(x))
#endif // APPLE
#endif // Windows

#endif // not MATLAB MEX
