#include <tbb/parallel_for.h>
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

    void parallel_for_c(void* kalman, void* indices, int length, size_t n, size_t block_size, void (*func)(void*, void*, int, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, block_size),
            [kalman, indices, length, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, indices, length, r.begin(), r.end());
            }
        );
    }
}
