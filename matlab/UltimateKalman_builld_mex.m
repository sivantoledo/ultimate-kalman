if ~exist('parallel','var')
    parallel = 'ultimate';
end

clear mex;
if ispc
    cd '..\c'
    if (~isempty(ver('MATLAB')))
        switch parallel
            case 'oddeven'
                disp('compiling and linking under MATLAB (parallel UltimateKalman) ...');
                mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c ultimatekalman_oddeven.c -lmwlapack -lmwblas
            case 'oddeven_nc'
                disp('compiling and linking under MATLAB (parallel UltimateKalman) ...');
                mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c ultimatekalman_oddeven_nc.c -lmwlapack -lmwblas
            case 'associative'
                disp('compiling and linking under MATLAB (parallel associative) ...');
                mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c kalman_associative_parallel.c -lmwlapack -lmwblas
            case 'ultimate'
                disp('compiling and linking under MATLAB ... (ultimate kalman)');
                mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c ultimatekalman.c flexible_arrays.c kalman_matrix_ops.c -lmwlapack -lmwblas
            case 'filter_smoother'
                disp('compiling and linking under MATLAB ... (kalman filter smoother)');
                mex -DNDEBUG -DBUILD_MEX -DBUILD_MATLAB -DBUILD_LAPACK_H -DBUILD_BLAS_H -DBUILD_WIN32_GETTIMEOFDAY ultimatekalmanmex.c kalman_filter_smoother.c -lmwlapack -lmwblas
            otherwise
                error(['unknown parallel variant ' parallel])
        end
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
    addpath '..\c'
end

if isunix
    cd '../c'
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
    addpath '../c'
end

if ismac
    warning('MacOS mex generation not yet tested');
    cd '../c'
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
    addpath '../c'
end


