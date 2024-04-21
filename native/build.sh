#!/bin/bash

gcc \
    -DBUILD_BLAS_UNDERSCORE -DBUILD_LAPACK_UNDERSCORE \
    -o rotation \
    rotation.c \
    ultimatekalman.c \
    -llapack -lblas -lm 

