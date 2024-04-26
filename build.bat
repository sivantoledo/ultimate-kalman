setlocal

echo COMPLILING JAVA

:: javac -cp commons-math3-3.6.1.jar -d bin -sourcepath src --target 1.8 --source 1.8 src\sivantoledo\kalman\*.java
javac -cp commons-math3-3.6.1.jar -d bin -sourcepath src --release 8 src\sivantoledo\kalman\*.java

echo CREATING JAR

echo Manifest-Version: 1.1> MANIFEST.TXT
echo Name: ultimatekalman>> MANIFEST.TXT
echo Implementation-Version: ultimatekalman jar produced %DATE% %TIME%>> MANIFEST.TXT
::echo Class-Path: ../jars/cmath.jar jars/log4j-api-2.2.jar ../jars/log4j-core-2.2.jar ../jars/opencsv-3.7.jar>> MANIFEST.TXT
del ultimatekalman.jar
cd bin
jar cvfm ..\ultimatekalman.jar ..\MANIFEST.TXT sivantoledo\* 
cd ..

java -cp ultimatekalman.jar;commons-math3-3.6.1.jar;bin sivantoledo.kalman.tests.Rotation 






