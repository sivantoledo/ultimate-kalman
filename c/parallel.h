#ifndef PARALLEL_H
#define PARALLEL_H

void parallel_set_thread_limit(int number_of_threads);
void parallel_set_blocksize   (int blocksize_in);

void foreach_in_range    (void (*func)(void*,        int, size_t, size_t), void* array ,               int length, size_t n);
void foreach_in_range_two(void (*func)(void*, void*, int, size_t, size_t), void* array1, void* array2, int length, size_t n);

void parallel_for_c_oddeven    (void* kalman, void* indices, int length, int** helper, size_t n, size_t block_size, void (*func)(void*, void*, int, int**, size_t, size_t));
void parallel_for_c_oddeven_nc (void* kalman, void* indices, int length,               size_t n, size_t block_size, void (*func)(void*, void*, int,        size_t, size_t));
void parallel_for_c_associative(void* kalman, void** helper, size_t l,                 size_t n, size_t block_size, void (*func)(void*, void**, size_t,    size_t, size_t));

void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride);

#endif