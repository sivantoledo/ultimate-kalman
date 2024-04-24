#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "provide an argument, the name of the C program you want to compile and run"
    exit
fi

echo generating test program $1

gcc \
    -O2 \
    -DBUILD_BLAS_UNDERSCORE -DBUILD_LAPACK_UNDERSCORE \
    -DBUILD_DEBUG_PRINTOUTSx \
    -o $1 \
    $1.c \
    ultimatekalman.c \
    -llapack -lblas -lm

echo generated test program

./$1

echo done running test

echo build script done


