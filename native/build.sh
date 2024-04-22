#!/bin/bash

gcc \
    -DBUILD_BLAS_UNDERSCORE -DBUILD_LAPACK_UNDERSCORE \
    -DBUILD_DEBUG_PRINTOUTSx \
    -o rotation \
    rotation.c \
    ultimatekalman.c \
    -llapack -lblas -lm 

