#include <cstddef>
#include <cstdint>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_scan.h>

// Define the parallel function with C linkage
extern "C" {
	static int nthreads  = 0;
	static int blocksize = 10;

	void parallel_set_thread_limit(int number_of_threads) {
		if (number_of_threads > 0) {
			nthreads = number_of_threads;
		}
	}

	void parallel_set_blocksize(int blocksize_in) {
		if (blocksize_in > 0) {
			blocksize = blocksize_in;
		}
	}

    void foreach_in_range(void (*func)(void*, int, size_t, size_t), void* array, int length, size_t n) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

    	//printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, blocksize),
            [array, length, func](const tbb::blocked_range<size_t>& subrange) {
                func(array, length, subrange.begin(), subrange.end());
            }
        );
    }
    
    void foreach_in_range_two(void (*func)(void*, void*, int, size_t, size_t), void* array1, void* array2, int length, size_t n) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

    	//printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, blocksize),
            [array1, array2, length, func](const tbb::blocked_range<size_t>& subrange) {
                func(array1, array2, length, subrange.begin(), subrange.end());
            }
        );
    }

#ifdef OBSOLETE
	// this is from oddeven
    void parallel_for_c_oddeven(void* kalman, void* indices, int length, int** helper, size_t n, size_t block_size, void (*func)(void*, void*, int, int**, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

    	//printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, blocksize>0 ? blocksize : block_size),
            [kalman, indices, length, helper, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, indices, length, helper, r.begin(), r.end());
            }
        );
    }
    
    // this is from oddeven_nc
    void parallel_for_c_oddeven_nc(void* kalman, void* indices, int length, size_t n, size_t block_size, void (*func)(void*, void*, int, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, n, block_size),
            [kalman, indices, length, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, indices, length, r.begin(), r.end());
            }
        );
    }
    
    // this is from associative
    void parallel_for_c_associative(void* kalman, void** helper, size_t l, size_t n, size_t block_size, void (*func)(void*, void**, size_t, size_t, size_t)) {
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);
        //tbb::global_control control(tbb::global_control::max_allowed_parallelism, NUM_OF_THREADS);

    	tbb::parallel_for(tbb::blocked_range<size_t>(0, n, block_size),
            [kalman, helper, l, func](const tbb::blocked_range<size_t>& r) {
                func(kalman, helper, l, r.begin(), r.end());
            }
        );
    }
#endif 
    
    void parallel_scan_c(void** input, void** sums, void* create_array , void* (*f)(void*, void*, void*, int, int), int length, int stride){
    	tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

        tbb::parallel_scan(
            tbb::blocked_range<size_t>(0, length),
            (void*) NULL, /* starting value (identity elements) */
            // now define the scan operation
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
            // and now the combining operation
            [f, create_array](void* left, void* right) { return f(left, right, create_array, -1, 0); }
            // there is also a version with an explicit is_final flag
        );
    }
    
}
