#ifndef KALMAN_TBB_H
#define KALMAN_TBB_H

int kalman_parallel_init(int number_of_threads);
int kalman_parallel_blocksize(int blocksize_in);

void foreach_in_range    (void* array ,               int length, size_t n, void (*func)(void*,        int, size_t, size_t));
void foreach_in_range_two(void* array1, void* array2, int length, size_t n, void (*func)(void*, void*, int, size_t, size_t));

void parallel_for_c_oddeven    (void* kalman, void* indices, int length, int** helper, size_t n, size_t block_size, void (*func)(void*, void*, int, int**, size_t, size_t));
void parallel_for_c_oddeven_nc (void* kalman, void* indices, int length,               size_t n, size_t block_size, void (*func)(void*, void*, int,        size_t, size_t));
void parallel_for_c_associative(void* kalman, void** helper, size_t l,                 size_t n, size_t block_size, void (*func)(void*, void**, size_t,    size_t, size_t));

void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride);

#endif