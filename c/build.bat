setlocal EnableDelayedExpansion

@ECHO ON

:: SET INT_TYPES=-DKALMAN_STEP_INDEX_TYPE_INT32 -DFARRAY_INDEX_TYPE_INT32 -DPARALLEL_INDEX_TYPE_INT32
SET INT_TYPES=-DKALMAN_STEP_INDEX_TYPE_INT64 -DFARRAY_INDEX_TYPE_INT64 -DPARALLEL_INDEX_TYPE_INT64
                     
:: possible mkl setvars args: ilp64
REM call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

SET C_FLAGS=/Zc:__cplusplus /Zp8 /GR /W3 /EHs /nologo /MD /O2 /Oy- -DNDEBUG

SET BLAS_LAPACK_FLAGS=-DBUILD_MKL

:: OpenMP (default) MKL, 64-bit pointers (lp64) but 32-bit row and column indices (ilp64 for 64-bit integers)
:: SET BLAS_LAPACK_LIBS=mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib
SET BLAS_LAPACK_LIBS=mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core_dll.lib libiomp5md.lib
:: Single-threaded library
SET BLAS_LAPACK_LIBS=mkl_intel_lp64.lib mkl_sequential.lib mkl_core_dll.lib 

SET ULTIMATE_C=kalman_ultimate.c ^
               kalman_conventional.c ^
               kalman_oddeven_smoother.c ^
               kalman_associative_smoother.c ^
               kalman_base.c ^
               kalman_explicit_representation.c ^
               matrix_ops.c ^
               flexible_arrays.c ^
               concurrent_set.c ^
               cmdline_args.c ^
               gettimeofday.c
               
SET CLIENTS=blastest rotation performance embarrassingly_parallel

ECHO UltimateKalman flags: %INT_TYPES%

del *.exe *.obj

set "ULTIMATE_O="
for %%F in (%ULTIMATE_C%) do (
  cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %INT_TYPES% %%F 
  
  set "name=%%F"
  set "name=!name:.c=.obj!"
  set "ULTIMATE_O=!ULTIMATE_O! !name!"
)

echo %ULTIMATE_O%

cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %INT_TYPES% parallel_sequential.c
cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %INT_TYPES% parallel_tbb.cpp 

for %%C in (%CLIENTS%) do (
    cl /c %C_FLAGS% -I. %BLAS_LAPACK_FLAGS% %INT_TYPES% %%C.c
	echo building %%C.exe
    cl %C_FLAGS% -Fe%%C.exe %ULTIMATE_O% parallel_sequential.obj %%C.obj %BLAS_LAPACK_LIBS% 
	echo building %%C_par.exe
    cl %C_FLAGS% -Fe%%C_par.exe %ULTIMATE_O% parallel_tbb.obj        %%C.obj %BLAS_LAPACK_LIBS% 
)

DEL *.obj
  
ECHO generated test programs

ECHO To run the generated binaries, invoke 
ECHO   "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
ECHO on the command line, to ensure that Windows can find the required DLLs.



