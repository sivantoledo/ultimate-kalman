# UltimateKalman: Flexible Kalman Filtering and Smoothing Using Orthogonal Transformations

This repository contains the source code of UltimateKalman in three different programming languages: MATLAB, C and Java.

Most of the documentation for UltimateKalman is available in an article [available on arXiv](https://arxiv.org/abs/2207.13526).

The code has been tested with MATLAB R2021b (version 9.11) and with GNU OCTAVE 7.1.0, both running under
Windows 11. 
MATLAB was configured to use Microsoft Visual C++ 2019 to compile C (mex) code. 
OCTAVE was configured to use mingw64 to compile C (mex) code. 

## Building and Testing the Code

Testing the MATLAB version is trivial. Launch MATLAB or OCTAVE, make sure that you are in the matlab directory
of this project (or that it is in MATLAB's search path), and run
`replication`.
MATLAB/OCTAVE will run a number
of tests and will produce the graphs in the article, except for the performance graphs. 

To build the Java version, run the Windows batch script 
`build.bat`.
The script invokes two command-like programs that are part of the Java Development Kit (JDK), `javac` and `jar`, 
so the JDK directory that contains these binaries must be on your path in order for the script to work.
The script is trivial and should be easy to port to Linux or MacOS.
The script builds an archive file called `ultimatekalman.jar`. To test the Java version, after
you have built it, run in MATLAB
`replication('Java')`.
It should produce exactly the same output as the MATLAB version.

Under OCTAVE, you may need to tell OCTAVE where to find the Java runtime by running in OCTAVE the command `setenv JAVA_HOME 'C:\Programs\jdk-11.0.2'` (replace the path with the path in which you installed the Java Runtime Environment or the Java Development Kit).

To build the C version, run the MATLAB script
`compile`.
It will compile the C version into a MATLAB-callable dynamic link library (a mex file).
To test the C version, run in MATLAB
`replication('C')`.
Again it should produce the same graphs.

The C code calls functions from two standard linear algebra libraries, the BLAS and LAPACK. These
libraries are standard, but their official interfaces are in Fortran (and for LAPACK, the implementation
of most of the functions is in fact in Fortran; the implementation of the BLAS is now often in C, 
but the interface we use is the Fotran interface). There is some variability in how these functions should
be called from C. For example, the matrix-multiplication routine may be `dgemm` or `dgemm_`, indexes
may be represented as 32-bit or 64-bit integers, and so on. To support this variability, the code
uses preprocessor variables that control how these functions are called. These variables are defined
appropriatedly in `compile.m` for the two environments in which I tested the code. If you move to another
platform, you may need to change them (e.g., BLAS and LAPACK functions have names that end with underscore under
MATLAB in Linux but not under MATLAB in Windows). 

These variables are:
- `BUILD_MEX` should be defined if building a mex file (MATLAB/OCTAVE callable library).
- `BUILD_MATLAB` should be defined if building a mex file in MATLAB.
- `BUILD_OCTAVE` should be defined if building a mex file in OCTAVE.
- `BUILD_BLAS_H` should be defined if the header `blas.h` is available.
- `BUILD_LAPACK_H` should be defined if the header `lapack.h` is available.
- `BUILD_BLAS_UNDERSCORE` should be defined if names of BLAS and LAPACK functions end with an underscore.
- `BUILD_BLAS_STRLEN_END` should be defined if the lengths of string arguments to BLAS and LAPACK
                              functions should be specified at the end of the argument list.
- `BUILD_WIN32_GETTIMEOFDAY` should be defined if the library function `gettimeofday` is not available
in the standard library.                    

Once you have built all three versions, you can produce the performance graphs shown in the 
article by running in MATLB
`replication('MATLAB',false,true)`.
**This will take a while!**
    
That's it!

To use the Java version with client code other than the MATLAB adapter class, simply include
`ultimatekalman.jar` and Apache Commons Math in the class path (this software comes with a particular
version of the Apache Commons Math library, `commons-math3-3.6.1.jar`).

To use the C version with client code other than the MATLAB adapter class, add to your project a single
C file, `ultimatekalman.c`, and a single header file, `ultimatekalman.h`.

## License

Copyright 2020-2022 Sivan Toledo.
 
 UltimateKalman is free software; you can redistribute it and/or modify
    it under the terms of either:

 the GNU Lesser General Public License as published by the Free
        Software Foundation; either version 3 of the License, or (at your
        option) any later version.

or

the GNU General Public License as published by the Free Software
        Foundation; either version 2 of the License, or (at your option) any
        later version.

or both in parallel, as here, 
    WITH THE ADDITIONAL REQUIREMENT 
    that if you use this software or derivatives of it, directly or indirectly, to produce
    research that is described in a research paper, you need to cite the most
    up-to-date version of the article that describes UltimateKalman in your paper.
    
Currently, the version to cite is [the version on arXiv](https://arxiv.org/abs/2207.13526).

UltimateKalman is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

You should have received copies of the GNU General Public License and the
    GNU Lesser General Public License along with this software.  If not,
    see https://www.gnu.org/licenses/.
    