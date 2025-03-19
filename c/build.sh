#!/bin/bash

# INT_TYPES="-DKALMAN_STEP_INDEX_TYPE_INT64 -DFARRAY_INDEX_TYPE_INT64 -DPARALLEL_INDEX_TYPE_INT64"
INT_TYPES="-DKALMAN_STEP_INDEX_TYPE_INT32 -DFARRAY_INDEX_TYPE_INT32 -DPARALLEL_INDEX_TYPE_INT32"

ARMPL_PATH="/opt/arm/armpl_24.10_gcc"
AMDPL_PATH="/specific/amd-gcc/5.0.0/gcc"
ONEAPI_PATH="/opt/intel/oneapi"

ULTIMATE_C="\
kalman_ultimate.c \
kalman_conventional.c \
kalman_oddeven_smoother.c \
kalman_associative_smoother.c \
kalman_base.c \
kalman_explicit_representation.c \
matrix_ops.c \
flexible_arrays.c \
concurrent_set.c \
cmdline_args.c"
ULTIMATE_O="${ULTIMATE_C//.c/.o}"

CLIENTS_C="blastest.c rotation.c performance.c embarrassingly_parallel.c"
CLIENTS_O="${CLIENTS_C//.c/.o}"
CLIENTS="${CLIENTS_C//.c/}"

case "$(uname)" in 
    Darwin)
        LIBDIR=""
        INCDIR="-I$(brew --prefix tbb)/include -DBUILD_BLAS_UNDERSCORE"
        SEQLIBS="-L$(brew --prefix openblas)/lib -llapack -lblas -lm"
        PARLIBS="-ltbbmalloc_proxy -ltbb -framework Accelerate -llapack -lblas -lm -lc"
        PARLIBS="-L$(brew --prefix openblas)/lib -llapack -lblas -lm -L$(brew --prefix tbb)/lib -ltbb -ltbbmalloc -ltbbmalloc_proxy"

        LIBDIR=""
        INCDIR="-I$(brew --prefix tbb)/include -DBUILD_BLAS_UNDERSCORE"
        SEQLIBS="-framework Accelerate -llapack -lblas -lm"
        PARLIBS="-framework Accelerate -llapack -lblas -lm -L$(brew --prefix tbb)/lib -ltbb -ltbbmalloc -ltbbmalloc_proxy"

        PRNLIBS=$PARLIBS
        ;;
    Linux)
	case "$(uname -m)" in
	    aarch64)
			echo "ARM"
			LIBDIR="-L${ARMPL_PATH}/lib/"
			INCDIR="-DBUILD_BLAS_UNDERSCORE"
			# sequential libraries; not sure why -lpthread was used, but it was included
			SEQLIBS="                  -larmpl_lp64 -lpthread -lm -ldl"
			PARLIBS="-ltbbmalloc_proxy -larmpl_lp64 -lpthread -lm -ldl -ltbb -lpthread -lm -ldl -ltbb"
			# parallel with nested TBB parallelism
			PRNLIBS="-ltbbmalloc_proxy -larmpl_lp64 -lpthread -lm -ldl -ltbb -lpthread -lm -ldl -ltbb"
	    ;;
	    x86_64)
		if grep -q "EPYC" /proc/cpuinfo; then
		   echo "Building for AMD EPYC"
		   LIBDIR="-L${AMDPL_PATH}/lib_LP64"
		   INCDIR="-DBUILD_BLAS_UNDERSCORE"
		   # sequential libraries; not sure why -lpthread was used, but it was included
		   SEQLIBS="                  -lflame -lblis-mt -laoclutils -lgomp -lpthread -lm -ldl"
		   PARLIBS="-ltbbmalloc_proxy -lflame -lblis-mt -laoclutils -lgomp -lpthread -lm -ldl -ltbb"
		   # parallel with nested TBB parallelism
		   PRNLIBS="-ltbbmalloc_proxy -lflame -lblis-mt -laoclutils -lpthread -lm -ldl -ltbb"
		else # assuming an Intel CPU
		    source ${ONEAPI_PATH}/setvars.sh
		    
		    LIBDIR=""
		    INCDIR="-DBUILD_MKL"
		    # sequential libraries; not sure why -lpthread was used, but it was included
		    SEQLIBS="                  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
		    PARLIBS="-ltbbmalloc_proxy -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -ltbb"
		    # parallel with nested TBB parallelism
		    PRNLIBS="-ltbbmalloc_proxy -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lpthread -lm -ldl -ltbb"
		fi
		;;
	    *)
		echo "I do not know how to build for this architecture (uname -m)"
		exit 1
		;;
       esac
       ;;
    *)
        echo "I do not know how to build the code on this operating system"
        exit 1
        ;;
esac

for C_SOURCE in $ULTIMATE_C; do
    echo compiling $C_SOURCE
    gcc -c -O2 $INCDIR $INT_TYPES $C_SOURCE
done

for C_SOURCE in $CLIENTS_C; do
    echo compiling $C_SOURCE
    gcc -c -O2 $INT_TYPES $C_SOURCE
done

echo compiling parallel_tbb.cpp
g++ -c -O2 $INCDIR $INT_TYPES -std=c++14 parallel_tbb.cpp

echo compiling parallel_sequential.c
gcc -c -O2 $INCDIR $INT_TYPES            parallel_sequential.c

echo LINKING

for CLIENT in $CLIENTS; do
    echo linking $CLIENT
    gcc $ULTIMATE_O parallel_sequential.o ${CLIENT}.o -o $CLIENT $LIBDIR $SEQLIBS
done

for CLIENT in $CLIENTS; do
    echo linking ${CLIENT}_par
    g++ $ULTIMATE_O parallel_tbb.o ${CLIENT}.o -o ${CLIENT}_par $LIBDIR $PARLIBS
done

case "$(uname)" in 
    Darwin)
        ;;
    Linux)
	case "$(uname -m)" in
	    aarch64)
			echo "Linux ARM, to set environment variables to run programs, run (under bash)"
			echo "export LD_LIBRARY_PATH=${ARMPL_PATH}/lib/"
			;;
	    x86_64)
			if grep -q "EPYC" /proc/cpuinfo; then
		 		echo "Linux AMD, to set environment variables to run programs, run (under bash)"
		 		echo "export LD_LIBRARY_PATH=${AMDPL_PATH}/lib_LP64"
  			else # assuming an Intel CPU
		  		echo "Linux Intel, to set environment variables to run programs, run (under bash)"
		  		echo "source ${ONEAPI_PATH}/setvars.sh"
			fi
			;;
	esac
	;;
esac

echo build script done


