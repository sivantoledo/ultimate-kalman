setlocal

 
@ECHO OFF
	
IF "%1"=="" (
  SET test=rotation
) ELSE (
  SET test=%1%
)

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
SET BLAS_LAPACK_LIBS=mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib
REM define MKL_ILP64 for 64-bit integers

REM -DBUILD_DEBUG_PRINTOUTS
REM -DBUILD_BLAS_STRLEN_END
REM -DBUILD_BLAS_UNDERSCORE

REM echo show %path% to see which directories store the MKL dll's
REM echo %path%

ECHO generating test program %test%.exe

:: kalman_ultimate.c ^
:: kalman_filter_smoother.c
:: kalman_oddeven.c

del %test%.exe

cl ^
  %C_FLAGS% ^
  -Fe%test%.exe ^
  -I. ^
  %BLAS_LAPACK_FLAGS% ^
  -DBUILD_WIN32_GETTIMEOFDAY ^
  kalman_ultimate.c ^
  kalman_filter_smoother.c ^
  kalman_parallel_sequential.c ^
  kalman_oddeven.c ^
  kalman_base.c ^
  matrix_ops.c ^
  flexible_arrays.c ^
  %test%.c ^
  %BLAS_LAPACK_LIBS% 

ECHO generated test program

DEL *.obj
  
%test%.exe

ECHO done running test 
  
ECHO build script done


