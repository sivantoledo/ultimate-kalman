setlocal

::set path=%path%;c:\Program Files\Java\jdk1.8.0_25\bin
::setenv JAVA_HOME c:\Program Files\Java\jdk1.8.0_25

echo COMPLILING JAVA

javac -cp commons-math3-3.6.1.jar -d bin -sourcepath src --target 1.8 --source 1.8 src\sivantoledo\kalman\*.java

echo CREATING JAR

echo Manifest-Version: 1.0> MANIFEST.TXT
echo Name: ultimatekalman>> MANIFEST.TXT
echo Implementation-Version: ultimatekalman jar produced %DATE% %TIME%>> MANIFEST.TXT
::echo Class-Path: ../jars/cmath.jar jars/log4j-api-2.2.jar ../jars/log4j-core-2.2.jar ../jars/opencsv-3.7.jar>> MANIFEST.TXT
del ultimatekalman.jar
cd bin
jar cvfm ..\ultimatekalman.jar ..\MANIFEST.TXT sivantoledo\* 
cd ..

del jss.zip
zip jss.zip src\sivantoledo\kalman\*.java commons-math3-3.6.1.jar gpl-3.0.txt lgpl-3.0.txt LICENSE.txt README.md native\*.c native\*.h matlab\*.m

goto :eof







