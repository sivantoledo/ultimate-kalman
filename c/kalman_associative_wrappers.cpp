#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
//#include "parallel_for_c.h"

// Define the parallel function with C linkage
extern "C" {
	static int nthreads = 0;

	int kalman_parallel_init(int number_of_threads) {
		if (number_of_threads > 0) {
			nthreads = number_of_threads;
		}
		return 0;
	}

	// not the same in in ultimatekalman_parallel ...
    void parallel_for_c(void* kalman, void** helper, size_t l, size_t n, size_t block_size, void (*func)(void*, void**, size_t, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);
        //tbb::global_control control(tbb::global_control::max_allowed_parallelism, NUM_OF_THREADS);

    	tbb::parallel_for(tbb::blocked_range<size_t>(0, n, block_size),
            [kalman, helper, l, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, helper, l, r.begin(), r.end());
            }
        );
    }

    /*
    void parallel_for_c(void* kalman, void* indices, int length, size_t n, size_t block_size, void (*func)(void*, void*, int, size_t, size_t)) {
    	if (nthreads>0) tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, block_size),
            [kalman, indices, length, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, indices, length, r.begin(), r.end());
            }
        );
    }
    */

    void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride){
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

        tbb::parallel_scan(
            tbb::blocked_range<size_t>(0, length),
            (void*) NULL,
            [input, sums, create_array, f, length, stride](const tbb::blocked_range<size_t>& r, void* sum, bool is_final_scan) {
                void* temp = sum;
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    int j = i + 1;
                    if (stride == -1){
                        j = length - 1 - i;
                    }

                    temp = f(temp, input[j], create_array, i, is_final_scan);

                    if (is_final_scan) {
                        sums[i] = temp;
                    }
                }
                return temp;
            },
            [f, create_array](void* left, void* right) { return f(left, right, create_array, -1, 0); }
        );
    }

}
