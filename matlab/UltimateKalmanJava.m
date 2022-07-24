classdef UltimateKalmanJava < handle
    % UltimateKalman   An implementation of the Paige-Saunders Kalman
    % filter and smoother by Sivan Toledo, Tel Aviv University.
    %
    % The filter is advanced by calling evolve and then observe every
    % step.
    %
    % To predict the next state(s) before providing the observations
    % (possibly before you have them) call observe and then filtered. Then
    % you can roll back and provide observations.
    %
    % UltimateKalman Methods:
    %    evolve   - Evolve the state using a linear matrix equation observe
    %    - Provide observations of the current state estimate - Return the
    %    most up to date estimate of a state vector forget   - Forget the
    %    oldest steps to save memory meory rollback - Roll back the object
    %    to an earlier step latest   - The index of the last step that was
    %    observed earliest - The index of the earliest step that has not
    %    been forgoten

    properties (Access = private)
        handle;
    end

    methods (Access = public)
        function kalman = UltimateKalmanJava()
            kalman = kalman@handle();
            kalman.handle = sivantoledo.kalman.UltimateKalman();
        end

        function delete(kalman)
            kalman.handle = -1; % invalid; not needed but added for safety
        end

        function i = earliest(kalman)
            %EARLIST   The index of the oldest step that has not been
            %          forgotten
            %
            %   kalman.EARLIEST() returns the index of the earliest step
            %   that has not been forgoten. Step numbers start from zero.
            i = kalman.handle.earliest();
        end

        function i = latest(kalman)
            %LATEST   The index of the last step that was observed.
            %
            %   kalman.LATEST() returns the index of the latest step that 
            %   was observed.
            i = kalman.handle.latest();
        end

        function evolve(kalman,n_i,H_i,F_i,c_i,K_i)
            %EVOLVE   Evolve the state using the given linear recurrence.
            %   kalman.EVOLVE(H_i,F_i,c_i,K_i) evolves the state using the recurrence
            %                H_i * u_i = F_i * u_{i-1} + c_i + epsilon
            %   where C_i is the covariance matrix of the error term epsilon.
            %   The covariance matrix C_i must be an instance of
            %   CovarianceMatrix. The matrices H and F are standard Matlab
            %   matrices, and be is a standard column vector.
            % 
            %   In the first step, this method does nothing; the first state is
            %   not an evolution of an earlier state. You can omit
            %   arguments or provide empty arguments in the first step.
            %
            %   The argument H can be empty, []. If it is, H is taken to be
            %   an identity matrix (possibly rectangular, if the dimension of
            %   this state is larger than the dimension of the previous state).
            if nargin==2 || isempty(F_i)
               kalman.handle.evolve(n_i);
               return;
            end
            if isempty(H_i)
                kalman.handle.evolve(n_i,          kalman.jmat(F_i),kalman.jvec(c_i),kalman.jcov(K_i));
            else
                kalman.handle.evolve(n_i,kalman.jmat(H_i),kalman.jmat(F_i),kalman.jvec(c_i),kalman.jcov(K_i));
            end
        end

        function observe(kalman,G_i,o_i,C_i)
            %OBSERVE   Provide observations of the current state.
            %   kalman.OBSERVE(G_i,o_i,C_i) provide observations that satisfy
            %   the linear equation
            %                o_i = G_i * u_i + delta_i
            %   where C_i is the covariance matrix of the error term delta_i.
            %   The covariance matrix C_i must be an instance of
            %   CovarianceMatrix. The matrix G_i is a standard Matlab
            %   matrix, and bo is a standard column vector. 
            % 
            %   kalman.OBSERVE() tells the algorithm that no
            %   observations are availble of the state of this step.
            % 
            %   This method must be called after advance and evolve.

            if nargin<=1
                kalman.handle.observe();
                return;
            end
            kalman.handle.observe(kalman.jmat(G_i),kalman.jvec(o_i),kalman.jcov(C_i));
        end

        %function [estimate,cov] = filtered(kalman)
            %FILTERED   Obtain an estimate of the current state and the 
            %           covariance of the estimate.
            %
            %   kalman.FILTERED() returns the estimate of the state of the
            %   last step that was observed (even if the observation was
            %   empty).
            %
            %   The method also returns the covariance of the estimate. The
            %   covariance is an instance of CovarianceMatrix.
            % 
            %   This method can only be called after observe.

        %    l = length(kalman.steps);
        %    estimate = kalman.steps{l}.Rdiag \ kalman.steps{l}.y;
        %    cov = kalman.steps{l}.estimatedCovariance;
        %end

        function [estimate,cov] = estimate(kalman,s)
            %ESTIMATE  Obtain the most recent estimate of the state of a step
            %          and the covariance of the estimate.
            %
            %   [estimate,cov] = kalman.ESTIMATE(s) returns an estimate of 
            %   the state of step s, which must still be in memory. It also 
            %   returns the covariance of the estimate.
            %
            %   [estimate,cov] = kalman.ESTIMATE() returns an estimate of
            %   the state of the latest step.
            %
            %   If kalman.smooth() was not called after step s was
            %   observed, the estimate is a filtered estimate. Otherwise it
            %   is the most recent smoothed estimate.
            
            %length(kalman.steps)
            if nargin<2
                s = -1;
            end
            jestimate = kalman.handle.estimate(s);
            if ~isempty(jestimate)
                estimate = jestimate.toArray();
                %class(estimate)
                if (~isnan(estimate(1)))
                  W        = kalman.handle.covarianceInverseFactor(s).getData();
                  cov = CovarianceMatrix(W,'W');
                else
                  cov = [];
                end
            else
                estimate = [];
                cov = [];
            end
        end

        function forget(kalman,s)
            %FORGET  Forget the oldest steps to save memory meory.
            %   
            %   kalman.FORGET(s) forgets steps s and older. 
            %
            %   kalman.FORGET()  forgets all but the last step.
            
            if nargin<2
                s = -1;
            end
            kalman.handle.forget(s);
        end

        function rollback(kalman,s)
            %ROLLBACK  Rolls the object back to just after the evolution
            %          to step s (but before the observation in step s).
            %
            %   kalman.ROLLBACK(s) rolls back the filter. Steps later than
            %   s are discareded, and so is the observation of step s, if
            %   there was any.
            % 
            %   This method allows you to roll back the object after it was
            %   used for prediction of future state that have not been
            %   observed yet. Typically steps s+1 and higher do not have
            %   any observation at the time of the rollback.
            
            if nargin<2
                s = -1;
            end
            kalman.handle.rollback(s);
        end
        
        function smooth(kalman)
            %SMOOTH  Compute smooth estimates of all the stored states.
            %
            %   kalman.SMOOTH() computes the smoothed estimated state of all
            %   the steps that are still in memory.
            % 
            %   This method must be called after the last step has been
            %   observed.
            kalman.handle.smooth();
        end

        function t = perftest(kalman,H,F,c,K,G,o,C,count,decimation)
            l = size(F,1);                                             % row dimension
            n = size(G,2);
            if size(H,1) ~= l                                          % this allows the user to pass [] for H
                if l == n
                    H = eye(l);
                else
                    H = [ eye(l) zeros(l,n - l)];
                end
            end

            [C_rep,C_type] = rep(C);
            [K_rep,K_type] = rep(K);
            %t = ultimatekalmanmex('perftest',kalman.handle,H,F,c,K_rep,double(K_type),G,o,C_rep,double(C_type),double(count),double(decimation));
        end


    end % methods

    methods (Static, Access = private)
        function jC = jcov(C)
            [R,type] = C.rep();
            switch type
                case 'w'
                    R = UltimateKalmanJava.jvec(R);
                    t = javaMethod('valueOf','sivantoledo.kalman.DiagonalCovarianceMatrix$Representation','DIAGONAL_INVERSE_STANDARD_DEVIATIONS');
                    jC = sivantoledo.kalman.DiagonalCovarianceMatrix(R,t);
                case 'W'
                    R = UltimateKalmanJava.jmat(R);
                    t = javaMethod('valueOf','sivantoledo.kalman.RealCovarianceMatrix$Representation','INVERSE_FACTOR');
                    jC = sivantoledo.kalman.RealCovarianceMatrix(R,t);
                case 'U'
                    R = UltimateKalmanJava.jmat(R);
                    t = javaMethod('valueOf','sivantoledo.kalman.RealCovarianceMatrix$Representation','FACTOR');
                    jC = sivantoledo.kalman.RealCovarianceMatrix(R,t);
                otherwise
                    type
                    error('cannot process covariance matrix type');
            end
        end

        function jA = jmat(A)
            jA = org.apache.commons.math3.linear.MatrixUtils.createRealMatrix(A);
        end

        function jv = jvec(v)
            jv = org.apache.commons.math3.linear.MatrixUtils.createRealVector(v);
        end

    end % static methods
end