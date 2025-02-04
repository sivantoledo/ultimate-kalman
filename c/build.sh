#!/bin/bash


if [ "$#" -ne 1 ]; then
    test=rotation
else
    test=$1
fi

echo generating test program $test

gcc \
    -O2 \
    -DBUILD_BLAS_UNDERSCORE -DBUILD_LAPACK_UNDERSCORE \
    -DBUILD_DEBUG_PRINTOUTSx \
    -o $test \
    $test.c \
    ultimatekalman.c \
    -framework Acceleratex \
    -llapack -lblas -lm

echo generated test program

./$test

echo done running test

echo build script done


