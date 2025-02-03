#!/bin/bash


if [ "$#" -ne 1 ]; then
    test=rotation
else
    test=$1
fi

echo run \"source /opt/intel/oneapi/setvars.sh\" under bash to set environment variables

# gcc -O2 -DBUILD_MKL -DBUILD_DEBUG_PRINTOUTSx -c performance.c ultimatekalman.c
gcc -O2 -c performance.c 

gcc -O2 -DBUILD_MKL -DNO_COVARIANCE_ESTIMATES \
                               -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_nc.o                  ultimatekalman.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman.o                     ultimatekalman.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_filter_smoother.o             kalman_filter_smoother.c

gcc -O2 -DBUILD_MKL -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven.o             ultimatekalman_oddeven.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_seq.o         ultimatekalman_oddeven.c

gcc -O2 -DBUILD_MKL -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_nc.o          ultimatekalman_oddeven_nc.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_nc_seq.o      ultimatekalman_oddeven_nc.c
	
gcc -O2 -DBUILD_MKL -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_associative.o                 kalman_associative.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_associative_seq.o             kalman_associative.c

gcc -O2 -DBUILD_MKL -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -DNDEBUG -c -o embarrassingly_parallel.o            embarrassingly_parallel.c
gcc -O2 -DBUILD_MKL            -DBUILD_DEBUG_PRINTOUTSx -DNDEBUG -c -o embarrassingly_parallel_seq.o        embarrassingly_parallel.c

g++ -std=c++11 -c ultimatekalman_oddeven_nc_wrappers.cpp
g++ -std=c++11 -c ultimatekalman_oddeven_wrappers.cpp
g++ -std=c++11 -c kalman_associative_wrappers.cpp

# -lmkl_tbb_thread
# -lmkl_sequential

g++ \
  -std=c++11 \
  -o embarrassingly_parallel \
  embarrassingly_parallel.o \
  ultimatekalman_oddeven_nc_wrappers.o \
  -ltbbmalloc_proxy \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl \
  -ltbb

g++ \
  -std=c++11 \
  -O2 -DBUILD_MKL \
  -o embarrassingly_parallel_seq \
  embarrassingly_parallel_seq.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lm -ldl

g++ \
  -std=c++11 \
  -o performance_associative \
  performance.o \
  kalman_associative.o \
  kalman_associative_wrappers.o \
  -ltbbmalloc_proxy \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl \
  -ltbb

g++ \
  -std=c++11 \
  -o performance_associative_seq \
  performance.o \
  kalman_associative_seq.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl

g++ \
  -std=c++11 \
  -o performance_oddeven_nc \
  performance.o \
  ultimatekalman_oddeven_nc.o \
  ultimatekalman_oddeven_nc_wrappers.o \
  -ltbbmalloc_proxy \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl \
  -ltbb

g++ \
  -std=c++11 \
  -o performance_oddeven_nc_seq \
  performance.o \
  ultimatekalman_oddeven_nc_seq.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl 

g++ \
  -std=c++11 \
  -o performance_oddeven \
  performance.o \
  ultimatekalman_oddeven.o \
  ultimatekalman_oddeven_wrappers.o \
  -ltbbmalloc_proxy \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl \
  -ltbb

g++ \
  -std=c++11 \
  -o performance_oddeven_seq \
  performance.o \
  ultimatekalman_oddeven_seq.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl 
  
# ----------- for testing nesting -------------------
g++ \
    -std=c++11 \
    -o performance_oddeven_nested \
    performance.o \
    ultimatekalman_oddeven.o \
    ultimatekalman_oddeven_wrappers.o \
    -ltbbmalloc_proxy \
    -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core \
    -lpthread \
    -lm -ldl \
    -ltbb

g++ \
    -std=c++11 \
    -o performance_oddeven_seq_nested \
    performance.o \
    ultimatekalman_oddeven_seq.o \
    -ltbbmalloc_proxy \
    -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core \
    -lpthread \
    -lm -ldl \
    -ltbb

# ---------------------------------------------------

g++ \
  -std=c++11 \
  -o performance_ultimate_nc \
  performance.o \
  ultimatekalman_nc.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl 

g++ \
  -std=c++11 \
  -o performance_ultimate \
  performance.o \
  ultimatekalman.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl 
  
g++ \
  -std=c++11 \
  -o performance_filter_smoother \
  performance.o \
  kalman_filter_smoother.o \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread \
  -lm -ldl 

#g++ \
#    -std=c++11 \
##    -o performance_parallelmkl \
#    performance.o \
#    ultimatekalman.o \
#    -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core \
#    -lpthread \
#    -lm -ldl \
#    -ltbb

#rm *.o

echo build script done


