#ifndef MEMORY_H
#define MEMORY_H

#ifdef PARALLEL

#ifdef MACOS

void* local_aligned_alloc(size_t alignment, size_t size) {
    void* p;
    int result = posix_memalign(&p, alignment, size);
    if (result != 0) return NULL;
    return p;
}
#define malloc(x) local_aligned_alloc(64,(x))

#else /* MACOS */

#define malloc(x) aligned_alloc(64,(x))

#endif /* not MACOS */

#endif /* PARALLEL */

#endif
