setlocal

set path=%path%;c:\Program Files\Java\jdk1.8.0_25\bin
setenv JAVA_HOME c:\Program Files\Java\jdk1.8.0_25

echo CREATING JAR

echo Manifest-Version: 1.0> MANIFEST.TXT
echo Name: ultimatekalman>> MANIFEST.TXT
echo Implementation-Version: ultimatekalman jar produced %DATE% %TIME%>> MANIFEST.TXT
::echo Class-Path: ../jars/cmath.jar jars/log4j-api-2.2.jar ../jars/log4j-core-2.2.jar ../jars/opencsv-3.7.jar>> MANIFEST.TXT
del ultimatekalman.jar
cd bin
jar cvfm ..\ultimatekalman.jar ..\MANIFEST.TXT sivantoledo\* 
cd ..

goto :eof







