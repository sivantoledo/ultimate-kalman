classdef KalmanPython < handle
    properties (Access = private)
        handle;
    end

    methods (Access = public)
        function kalman = KalmanPython(options)
            KalmanPython.initializeClass();

            if nargin < 1 || ~isfield(options,'algorithm')
                error('options must includes at least an algorithms field');
            end
            kalman = kalman@handle();
            disp('calling factory')
            disp(options.algorithm)
            factory = py.kalman_factory.kalman_factory(options.algorithm, options);
            kalman.handle = factory();
        end

        function delete(kalman)
            kalman.handle = [];
        end

        function i = earliest(kalman)
            i = double(kalman.handle.earliest());
        end

        function i = latest(kalman)
            i = double(kalman.handle.latest());
        end

        function evolve(kalman, n_i, H_i, F_i, c_i, K_i)
            n_i = int64(n_i);
            if nargin == 2 || isempty(F_i)
                kalman.handle.evolve(n_i,[],[],[],[]);
                return;
            end
            if isempty(H_i)
                H_i = eye(n_i);
            end
            kalman.handle.evolve(n_i, kalman.pymat(H_i), kalman.pymat(F_i), kalman.pyvec(c_i), kalman.pycov(K_i));
        end

        function observe(kalman, G_i, o_i, C_i)
            if nargin <= 1
                kalman.handle.observe();
                return;
            end
            kalman.handle.observe(kalman.pymat(G_i), kalman.pyvec(o_i), kalman.pycov(C_i));
        end

        function [estimate, cov] = estimate(kalman, s)
            if nargin < 2
                jestimate = kalman.handle.estimate();
            else
                jestimate = kalman.handle.estimate(s);
            end
            rep = jestimate{2}.rep();
            cov     = double(rep{1});
            covType = rep{2};
            if ~isempty(jestimate)
                estimate = double(jestimate{1})';
                cov = CovarianceMatrix(cov,covType);
            else
                estimate = [];
                cov = [];
            end
        end

        function forget(kalman, s)
            if nargin < 2
                s = -1;
            end
            fprintf('forget %d %d\n',s,kalman.latest())
            kalman.handle.forget(int64(s));
        end

        function rollback(kalman, s)
            if nargin < 2
                s = -1;
            end
            kalman.handle.rollback(int64(s));
        end

        function smooth(kalman)
            kalman.handle.smooth();
        end
    end

    methods (Static, Access = private)
        function initializeClass()
            persistent isInitialized
            if ~isempty(isInitialized)
                return
            end

            fprintf('performing one-time initialization of KalmanPython\n')
            currentFile = mfilename('fullpath');
            [currentDir, ~, ~] = fileparts(currentFile);
            pythondir = fullfile(currentDir, '../python');
            py.sys.path().append(pythondir);
            isInitialized = true;
        end

        function pyC = pycov(C)
            [cov,type] = C.rep();
            pyCov = py.numpy.reshape(py.numpy.array(cov), int32([size(cov,1), size(cov,2)]));
            pyC = py.covariance_matrix.CovarianceMatrix(pyCov,type);
        end

        function pyA = pymat(A)
            %pyA = py.numpy.array(A);
            pyA = py.numpy.reshape(py.numpy.array(A), int32([size(A,1), size(A,2)]));
        end

        function pyv = pyvec(v)
            %pyv = py.numpy.array(v(:));
            pyv = py.numpy.reshape(py.numpy.array(v), int32([size(v,1), size(v,2)]));
        end
    end
end
