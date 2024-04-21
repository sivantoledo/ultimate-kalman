setlocal

@ECHO ON

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

::echo %MKLROOT%
::echo %INCLUDE%

:: possible mkl setvars args: ilp64
 
@ECHO ON
	
IF "%1"=="" (
  ECHO no arguments, building everything
  SET dll=true
  SET test=true
) 

IF "%1"=="dll"       (SET dll=true) 
IF "%1"=="test"      (SET test=true) 

:: /Og rather than /Ox ?

IF "%dll%"=="true" (
  ECHO generating library dll
  cl ^
  /O2 /Ox /EHsc ^
  -D_JNI_IMPLEMENTATION_ ^
  -DINPROCESS ^
  -I "%JAVA_HOME%\include" ^
  -I "%JAVA_HOME%\include\win32" ^
  -Fejni-sdcard.dll ^
    sdcard\tau_atlas_logger_SDCardController.c  ^
    /MD /LD
  
  COPY jni-sdcard.dll %ATLAS_DIST%\windows
  ECHO generated sdcard dll
)


IF "%test%"=="true" (
  ECHO generating test program
  cl ^
  -Ferotation.exe ^
  -I. ^
  -DBUILD_WIN32_GETTIMEOFDAY ^
  -DBUILD_MKL_H ^
  rotation.c ^
  ultimatekalman.c ^
  mkl_rt.lib ^
  /MT



  
  ECHO generated test
)


ECHO done


