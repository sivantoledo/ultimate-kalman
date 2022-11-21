clear mex;
cd '..\native'
if (~isempty(ver('MATLAB')))
  disp('compiling and linking under MATLAB ...');
  mex -DBUILD_MEX ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
  %mex -DBUILD_MEX ultimatekalmanmex.c ultimatekalman.c 
end
if (~isempty(ver('Octave')))
  disp('compiling and linking under Octave ...');
  mkoctfile --mex -v -fmax-errors=5 '-D BUILD_MEX' '-D BUILD_OCTAVE' '-I/Programs/OpenBLAS-0.3.21-x86/include' ultimatekalmanmex.c ultimatekalman.c -llibopenblas -lliblapack
  %mex -DBUILD_MEX ultimatekalmanmex.c ultimatekalman.c 
end
disp('compiling and linking done; now testing');
cd '..\matlab'
addpath '..\native'
