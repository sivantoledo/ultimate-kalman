#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
//#include "parallel_for_c.h"

// Define the parallel function with C linkage
extern "C" {
	static int nthreads  = 0;
	static int blocksize = 0;

	int kalman_parallel_init(int number_of_threads) {
		if (number_of_threads > 0) {
			nthreads = number_of_threads;
		}
		return 0;
	}

	int kalman_parallel_blocksize(int blocksize_in) {
		if (blocksize_in > 0) {
			blocksize = blocksize_in;
		}
		return 0;
	}


    void parallel_for_c(void* kalman, void* indices, int length, int** helper, size_t n, size_t block_size, void (*func)(void*, void*, int, int**, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

    	//printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, blocksize>0 ? blocksize : block_size),
            [kalman, indices, length, helper, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, indices, length, helper, r.begin(), r.end());
            }
        );
    }
}
