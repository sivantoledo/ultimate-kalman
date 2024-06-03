clear mex;
if ispc
    cd '..\native'
    if (~isempty(ver('MATLAB')))
        disp('compiling and linking under MATLAB ...');
        mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
    end
    if (~isempty(ver('Octave')))
        disp('compiling and linking under Octave ...');
        mkoctfile --mex -v -fmax-errors=5 ...
            -DNDEBUG ...
            -DBUILD_MEX ...
            -DBUILD_OCTAVE ...
            -DBUILD_BLAS_UNDERSCORE ...
            -DBUILD_BLAS_STRLEN_END ...
            ultimatekalmanmex.c ultimatekalman.c -llibopenblas -lliblapack
    end
    disp('compiling and linking done');
    cd '..\matlab'
    addpath '..\native'
end

if isunix
    cd '../native'
    if (~isempty(ver('MATLAB')))
        disp('compiling and linking under MATLAB ...');
        mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
    end
    if (~isempty(ver('Octave')))
        disp('compiling and linking under Octave ...');
        mkoctfile --mex -v -fmax-errors=5 ...
            -DNDEBUG ...
            -DBUILD_MEX ...
            -DBUILD_OCTAVE ...
            -DBUILD_BLAS_UNDERSCORE ...
            -DBUILD_BLAS_STRLEN_END ...
            ultimatekalmanmex.c ultimatekalman.c -llibopenblas -lliblapack
    end
    disp('compiling and linking done');
    cd '../matlab'
    addpath '../native'
end

if ismac
    warning('MacOS mex generation not yet tested');
    cd '../native'
    if (~isempty(ver('MATLAB')))
        disp('compiling and linking under MATLAB ...');
        mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H ultimatekalmanmex.c ultimatekalman.c -lmwlapack -lmwblas
    end
    if (~isempty(ver('Octave')))
        disp('compiling and linking under Octave ...');
        mkoctfile --mex -v -fmax-errors=5 ...
            -DNDEBUG ...
            -DBUILD_MEX ...
            -DBUILD_OCTAVE ...
            -DBUILD_BLAS_UNDERSCORE ...
            -DBUILD_BLAS_STRLEN_END ...
            ultimatekalmanmex.c ultimatekalman.c -llibopenblas -lliblapack
    end
    disp('compiling and linking done');
    cd '../matlab'
    addpath '../native'
end


