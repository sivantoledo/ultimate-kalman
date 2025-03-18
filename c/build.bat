setlocal

@ECHO ON

REM call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

::echo %MKLROOT%
::echo %INCLUDE%

:: possible mkl setvars args: ilp64

:: /Og rather than /Ox ?

SET C_FLAGS=/Zc:__cplusplus /Zp8 /GR /W3 /EHs /nologo /MD /O2 /Oy- -DNDEBUG

SET BLAS_LAPACK_FLAGS=-DBUILD_MKL

:: OpenMP (default) MKL
:: SET BLAS_LAPACK_LIBS=mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib
SET BLAS_LAPACK_LIBS=mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core_dll.lib libiomp5md.lib
:: Single-threaded library
SET BLAS_LAPACK_LIBS=mkl_intel_lp64.lib mkl_sequential.lib mkl_core_dll.lib 
REM define MKL_ILP64 for 64-bit integers

::SET ULTIMATE_FLAGS=-DBUILD_WIN32_GETTIMEOFDAY 
SET ULTIMATE_FLAGS=-DKALMAN_STEP_INDEX_TYPE_INT32 -DFARRAY_INDEX_TYPE_INT32 -DPARALLEL_INDEX_TYPE_INT32
SET ULTIMATE_FLAGS=-DKALMAN_STEP_INDEX_TYPE_INT64 -DFARRAY_INDEX_TYPE_INT64 -DPARALLEL_INDEX_TYPE_INT64
                     
SET ULTIMATE_OBJECTS=kalman_ultimate.obj ^
                     kalman_conventional.obj ^
                     kalman_oddeven_smoother.obj ^
                     kalman_associative_smoother.obj ^
                     kalman_base.obj ^
                     kalman_explicit_representation.obj ^
                     matrix_ops.obj ^
                     flexible_arrays.obj ^
                     concurrent_set.obj ^
                     cmdline_args.obj ^
                     gettimeofday.obj
                     
REM -DBUILD_DEBUG_PRINTOUTS
REM -DBUILD_BLAS_STRLEN_END
REM -DBUILD_BLAS_UNDERSCORE

REM echo show %path% to see which directories store the MKL dll's
REM echo %path%

ECHO generating test program %test%.exe

ECHO UltimateKalman flags: %ULTIMATE_FLAGS%

:: kalman_ultimate.c ^
:: kalman_filter_smoother.c
:: kalman_oddeven.c

:: del %test%.exe
del *.exe *.obj

cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_ultimate.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_conventional.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_oddeven_smoother.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_associative_smoother.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_base.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% kalman_explicit_representation.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% matrix_ops.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% flexible_arrays.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% concurrent_set.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% parallel_sequential.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% gettimeofday.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% cmdline_args.c 

cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% parallel_tbb.cpp 

cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% embarrassingly_parallel.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% performance.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% rotation.c 
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %ULTIMATE_FLAGS% blastest.c 

cl %C_FLAGS% -Feperformance.exe             %ULTIMATE_OBJECTS% parallel_sequential.obj performance.obj            %BLAS_LAPACK_LIBS% 
cl %C_FLAGS% -Ferotation.exe                %ULTIMATE_OBJECTS% parallel_sequential.obj rotation.obj               %BLAS_LAPACK_LIBS% 
cl %C_FLAGS% -Feblastest.exe                %ULTIMATE_OBJECTS% parallel_sequential.obj blastest.obj               %BLAS_LAPACK_LIBS% 

cl %C_FLAGS% -Feperformance_par.exe         %ULTIMATE_OBJECTS% parallel_tbb.obj        performance.obj             %BLAS_LAPACK_LIBS% 
cl %C_FLAGS% -Feembarrassingly_parallel.exe %ULTIMATE_OBJECTS% parallel_tbb.obj        embarrassingly_parallel.obj %BLAS_LAPACK_LIBS% 

DEL *.obj
  
ECHO generated test programs

ECHO To run the generated binaries, invoke 
ECHO   "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
ECHO on the command line, to ensure that Windows can find the required DLLs.



