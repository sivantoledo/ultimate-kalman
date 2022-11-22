clear mex;
cd '..\native'
if (~isempty(ver('MATLAB')))
  disp('compiling and linking under MATLAB ...');
  mex -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
end
if (~isempty(ver('Octave')))
  disp('compiling and linking under Octave ...');
  mkoctfile --mex -v -fmax-errors=5 ...
    -DBUILD_MEX ...
    -DBUILD_OCTAVE ...
    -DBUILD_BLAS_UNDERSCORE ...
    -DBUILD_BLAS_STRLEN_END ...
    ultimatekalmanmex.c ultimatekalman.c -llibopenblas -lliblapack
end
disp('compiling and linking done; now testing');
cd '..\matlab'
addpath '..\native'
