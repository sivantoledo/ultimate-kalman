setlocal

::set ATLAS_DIST=C:\files\atlas\atlas-distribution
::set ATLAS_JAVA=C:\files\atlas\atlas-java
::set ATLAS_DEPS=C:\files\native\atlas-deps-native

::set path=c:\programs\jdk-11.0.2\bin;%path%
::set JAVA_HOME=c:\programs\jdk-11.0.2

::ECHO building header files

::javac -h dsp -d /tmp/junk -cp ^

@ECHO ON

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
 
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
  cl -I. ultimate-kalman.c /MT
  
::  COPY jni-sdcard.dll %ATLAS_DIST%\windows
  ECHO generated test
)


::cl ^
::  ...
::  -Feserver-usrp.exe ^
::    radio\server.c ^
::  /MT ^
::  /link ^
::  /LIBPATH:"c:\Program Files (x86)\UHD\lib" uhd.lib ^  
::)

:done

ECHO done


