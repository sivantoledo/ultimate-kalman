#!/bin/bash


if [ "$#" -ne 1 ]; then
    test=rotation
else
    test=$1
fi

echo generating test program $test

case "$(uname)" in 
    Darwin)
        LIBDIR="-framework Accelerate"
        ;;
    Linux)
        LIBDIR=""
        ;;
    *)
        echo "I do not know how to build the code on this operating system"
        exit 1
        ;;
esac


gcc \
    -O2 \
    -DBUILD_BLAS_UNDERSCORE -DBUILD_LAPACK_UNDERSCORE \
    -DBUILD_DEBUG_PRINTOUTSx \
    -o $test \
    $test.c \
    ultimatekalman.c \
    kalman_base.c \
    kalman_parallel_sequential.c \
    kalman_matrix_ops.c \
    flexible_arrays.c \
    $LIBDIR \
    -llapack -lblas -lm

echo generated test program

./$test

echo done running test

echo build script done


