/*
 * parallel_tbb.cpp
 *
 * A TBB-based implementation of the parallel primitives.
 *
 * Copyright (c) 2024-2025 Sivan Toledo and Shahaf Gargir
 */

#include <cstddef>
#include <cstdint>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_scan.h>
#include <tbb/spin_mutex.h>

extern "C" {

#include "parallel.h"
#include "concurrent_set.h"

struct spin_mutex_st {
  tbb::spin_mutex mutex;
};

static int nthreads = 0;
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

void foreach_in_range(void (*func)(void*, parallel_index_t, parallel_index_t, parallel_index_t), void* array, parallel_index_t length, parallel_index_t n) {
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

  //printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
    tbb::parallel_for(tbb::blocked_range<parallel_index_t>(0, n, blocksize),
        [array, length, func](const tbb::blocked_range<size_t>& subrange) {
          func(array, length, subrange.begin(), subrange.end());
        }
    );
  }

  void foreach_in_range_two(void (*func)(void*, void*, parallel_index_t, parallel_index_t, parallel_index_t), void* array1, void* array2, parallel_index_t length, parallel_index_t n) {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nthreads);

    //printf("blocksize = %d\n",blocksize>0 ? blocksize : block_size);
    tbb::parallel_for(tbb::blocked_range<parallel_index_t>(0, n, blocksize),
        [array1, array2, length, func](const tbb::blocked_range<size_t>& subrange) {
          func(array1, array2, length, subrange.begin(), subrange.end());
        }
    );
  }

  spin_mutex_t* spin_mutex_create() {
    spin_mutex_t* wrapper = (spin_mutex_t*) malloc(sizeof(spin_mutex_t));
    if (wrapper) {
      // Construct the spin_mutex in place
      new (&wrapper->mutex) tbb::spin_mutex();
    }
    return wrapper;
  }

  void spin_mutex_lock(spin_mutex_t* mutex) {
    if (mutex) {
      mutex->mutex.lock();
    }
  }

  void spin_mutex_unlock(spin_mutex_t* mutex) {
    if (mutex) {
      mutex->mutex.unlock();
    }
  }

  void spin_mutex_destroy(spin_mutex_t* mutex) {
    if (mutex) {
      mutex->mutex.~spin_mutex();
      free(mutex);
    }
  }

}

/******************************************************************************/
/* END OF FILE                                                                */
/******************************************************************************/
