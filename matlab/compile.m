clear mex;
cd '..\native'
disp('compiling and linking...');
mex -DBUILD_MEX ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
%mex -DBUILD_MEX ultimatekalmanmex.c ultimatekalman.c 
disp('compiling and linking done; now testing');
cd '..\matlab'
addpath 'C:\Users\stoledo\git\ultimate-kalman\native'
