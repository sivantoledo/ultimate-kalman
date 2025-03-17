#!/bin/bash


if [ "$#" -ne 1 ]; then
    test=rotation
else
    test=$1
fi

# On Intel servers, we installed Intel's oneAPI package, which includes both TBB and MKL (which includes optimized BLAS and LAPACK).
# On ARM servers (Graviton3), we installed tbb from Ubuntu's libtbb-dev package, and the ARM Performance Libraries (downloaded as a tar file).

echo On Linux Intel, run \"source /opt/intel/oneapi/setvars.sh\" under bash to set environment variables
echo On Linux ARM, run \"export LD_LIBRARY_PATH=/opt/arm/armpl_24.10_gcc/lib/\"

ULTIMATE_C="\
kalman_ultimate.c \
kalman_conventional.c \
kalman_oddeven_smoother.c \
kalman_associative_smoother.c \
kalman_base.c \
kalman_explicit_representation.c \
matrix_ops.c \
flexible_arrays.c \
concurrent_set.c"

CLIENTS_C="blastest.c rotation.c performance.c"

ULTIMATE_O="${ULTIMATE_C//.c/.o}"
CLIENTS_O="${CLIENTS_C//.c/.o}"
CLIENTS="${CLIENTS_C//.c/}"

INT_TYPES="-DKALMAN_STEP_INDEX_TYPE_INT64 -DFARRAY_INDEX_TYPE_INT64 -DPARALLEL_INDEX_TYPE_INT64"
INT_TYPES="-DKALMAN_STEP_INDEX_TYPE_INT32 -DFARRAY_INDEX_TYPE_INT32 -DPARALLEL_INDEX_TYPE_INT32"

case "$(uname)" in 
    Darwin)
        LIBDIR=""
        INCDIR="-I$(brew --prefix tbb)/include -DBUILD_BLAS_UNDERSCORE -DMACOS"
        SEQLIBS="-L$(brew --prefix openblas)/lib -llapack -lblas -lm"
        PARLIBS="-ltbbmalloc_proxy -ltbb -framework Accelerate -llapack -lblas -lm -lc"
        PARLIBS="-L$(brew --prefix openblas)/lib -llapack -lblas -lm -L$(brew --prefix tbb)/lib -ltbb -ltbbmalloc -ltbbmalloc_proxy"

        LIBDIR=""
        INCDIR="-I$(brew --prefix tbb)/include -DBUILD_BLAS_UNDERSCORE -DMACOS"
        SEQLIBS="-framework Accelerate -llapack -lblas -lm"
        PARLIBS="-framework Accelerate -llapack -lblas -lm -L$(brew --prefix tbb)/lib -ltbb -ltbbmalloc -ltbbmalloc_proxy"

        PRNLIBS=$PARLIBS
        ;;
    Linux)
	case "$(uname -m)" in
	    aarch64)
		echo "ARM"
		LIBDIR="-L/opt/arm/armpl_24.10_gcc/lib/"
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
		   LIBDIR="-L/specific/amd-gcc/5.0.0/gcc/lib_LP64"
		   INCDIR="-DBUILD_BLAS_UNDERSCORE"
		   # sequential libraries; not sure why -lpthread was used, but it was included
		   SEQLIBS="                  -lflame -lblis-mt -laoclutils -lgomp -lpthread -lm -ldl"
		   PARLIBS="-ltbbmalloc_proxy -lflame -lblis-mt -laoclutils -lgomp -lpthread -lm -ldl -ltbb"
		   # parallel with nested TBB parallelism
		   PRNLIBS="-ltbbmalloc_proxy -lflame -lblis-mt -laoclutils -lpthread -lm -ldl -ltbb"
		else # assuming an Intel CPU
		    source /opt/intel/oneapi/setvars.sh
		    
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

exit


gcc -O2 $INT_TYPES -c blastest.c 
gcc -O2 $INT_TYPES -c performance.c 
gcc -O2 $INT_TYPES -c rotation.c 

gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_explicit_representation.o     kalman_explicit_representation.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_base.o                        kalman_base.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o matrix_ops.o                         matrix_ops.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o flexible_arrays.o                    flexible_arrays.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o concurrent_set.o                     concurrent_set.c

gcc -O2 $INCDIR $INT_TYPES  -DBUILD_DEBUG_PRINTOUTSx -c -o parallel_sequential.o                parallel_sequential.c
g++     $INCDIR $INT_TYPES  -std=c++11               -c                                         parallel_tbb.cpp

gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_ultimate.o                    kalman_ultimate.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_conventional.o                kalman_conventional.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_oddeven_smoother.o            kalman_oddeven_smoother.c
gcc -O2 $INCDIR $INT_TYPES -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_associative_smoother.o        kalman_associative_smoother.c



#gcc -O2 $INCDIR -DNO_COVARIANCE_ESTIMATES \
#                           -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_nc.o                  ultimatekalman.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman.o                     ultimatekalman.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_conventional.o             kalman_conventional.c

#gcc -O2 $INCDIR -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother.o             ultimatekalman_oddeven_smoother.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother_seq.o         ultimatekalman_oddeven_smoother.c

#gcc -O2 $INCDIR -DNO_COVARIANCE_ESTIMATES \
#                           -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_nc.o                  ultimatekalman.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman.o                     ultimatekalman.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_conventional.o             kalman_conventional.c

#gcc -O2 $INCDIR -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother.o             ultimatekalman_oddeven_smoother.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother_seq.o         ultimatekalman_oddeven_smoother.c

#gcc -O2 $INCDIR -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother_nc.o          ultimatekalman_oddeven_smoother_nc.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o ultimatekalman_oddeven_smoother_nc_seq.o      ultimatekalman_oddeven_smoother_nc.c
	
#gcc -O2 $INCDIR -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_associative_smoother.o                 kalman_associative_smoother.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -c -o kalman_associative_smoother_seq.o             kalman_associative_smoother.c

#gcc -O2 $INCDIR -DPARALLEL -DBUILD_DEBUG_PRINTOUTSx -DNDEBUG -c -o embarrassingly_parallel.o            embarrassingly_parallel.c
#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -DNDEBUG -c -o embarrassingly_parallel_seq.o        embarrassingly_parallel.c

#gcc -O2 $INCDIR            -DBUILD_DEBUG_PRINTOUTSx -DNDEBUG -c -o kalman_parallel_sequential.o         kalman_parallel_sequential.c

#g++ -std=c++11 $INCDIR -c ultimatekalman_oddeven_smoother_nc_wrappers.cpp
#g++ -std=c++11 $INCDIR -c ultimatekalman_oddeven_smoother_wrappers.cpp
#g++ -std=c++11 $INCDIR -c kalman_associative_smoother_wrappers.cpp

# -lmkl_tbb_thread
# -lmkl_sequential

echo LINKING

g++ \
    -std=c++11 \
    -o rotation_seq \
    rotation.o \
    kalman_ultimate.o \
    kalman_conventional.o \
    kalman_oddeven_smoother_seq.o \
    kalman_associative_smoother_seq.o \
    kalman_base.o \
    flexible_arrays.o \
    concurrent_set.o \
    matrix_ops.o \
    parallel_sequential.o \
    $LIBDIR $SEQLIBS

g++ \
    -std=c++11 \
    -o performance_seq \
    performance.o \
    kalman_ultimate.o \
    kalman_conventional.o \
    kalman_oddeven_smoother_seq.o \
    kalman_associative_smoother_seq.o \
    kalman_base.o \
    flexible_arrays.o \
    concurrent_set.o \
    matrix_ops.o \
    parallel_sequential.o \
    $LIBDIR $SEQLIBS

g++ \
    -std=c++11 \
    -o rotation_par \
    rotation.o \
    kalman_ultimate.o \
    kalman_conventional.o \
    kalman_oddeven_smoother.o \
    kalman_associative_smoother.o \
    kalman_base_par.o \
    flexible_arrays.o \
    concurrent_set.o \
    matrix_ops.o \
    parallel_tbb.o \
    $LIBDIR $PARLIBS


g++ \
    -std=c++11 \
    -o performance_par \
    performance.o \
    kalman_ultimate.o \
    kalman_conventional.o \
    kalman_oddeven_smoother.o \
    kalman_associative_smoother.o \
    kalman_base_par.o \
    flexible_arrays.o \
    concurrent_set.o \
    matrix_ops.o \
    parallel_tbb.o \
    $LIBDIR $PARLIBS

echo DONE

exit

g++ \
  -std=c++11 \
  -o embarrassingly_parallel \
  embarrassingly_parallel.o \
  ultimatekalman_oddeven_smoother_nc_wrappers.o \
  $LIBDIR $PARLIBS

g++ \
  -std=c++11 \
  -o embarrassingly_parallel_seq \
  embarrassingly_parallel_seq.o \
  $LIBDIR $SEQLIBS

g++ \
  -std=c++11 \
  -o performance_associative_smoother \
  performance.o \
  kalman_associative_smoother.o \
  kalman_associative_smoother_wrappers.o \
  $LIBDIR $PARLIBS

g++ \
  -std=c++11 \
  -o performance_associative_smoother_seq \
  performance.o \
  kalman_associative_smoother_seq.o \
  $LIBDIR $SEQLIBS

g++ \
  -std=c++11 \
  -o performance_oddeven_smoother_nc \
  performance.o \
  ultimatekalman_oddeven_smoother_nc.o \
  ultimatekalman_oddeven_smoother_nc_wrappers.o \
  $LIBDIR $PARLIBS

g++ \
  -std=c++11 \
  -o performance_oddeven_smoother_nc_seq \
  performance.o \
  ultimatekalman_oddeven_smoother_nc_seq.o \
  $LIBDIR $SEQLIBS

g++ \
  -std=c++11 \
  -o performance_oddeven_smoother \
  performance.o \
  ultimatekalman_oddeven_smoother.o \
  ultimatekalman_oddeven_smoother_wrappers.o \
  $LIBDIR $PARLIBS

g++ \
  -std=c++11 \
  -o performance_oddeven_smoother_seq \
  performance.o \
  ultimatekalman_oddeven_smoother_seq.o \
  $LIBDIR $SEQLIBS
  
# ----------- for testing nesting -------------------
g++ \
    -std=c++11 \
    -o performance_oddeven_smoother_nested \
    performance.o \
    ultimatekalman_oddeven_smoother.o \
    ultimatekalman_oddeven_smoother_wrappers.o \
    $LIBDIR $PRNLIBS

g++ \
    -std=c++11 \
    -o performance_oddeven_smoother_seq_nested \
    performance.o \
    ultimatekalman_oddeven_smoother_seq.o \
    $LIBDIR $PRNLIBS

# ---------------------------------------------------

g++ \
  -std=c++11 \
  -o performance_ultimate_nc \
  performance.o \
  ultimatekalman_nc.o \
  $LIBDIR $SEQLIBS

g++ \
  -std=c++11 \
  -o performance_ultimate \
  performance.o \
  ultimatekalman.o \
  $LIBDIR $SEQLIBS
  
g++ \
  -std=c++11 \
  -o performance_conventional \
  performance.o \
  kalman_conventional.o \
  $LIBDIR $SEQLIBS

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


