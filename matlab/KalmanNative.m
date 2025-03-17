classdef UltimateKalmanNative < handle
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
    %    evolve   - Evolve the state using a linear matrix equation
    %    observe  - Provide observations of the current state
    %    estimate - Return the most up to date estimate of a state vector
    %    forget   - Forget the oldest steps to save memory meory
    %    rollback - Roll back the object to an earlier step
    %    latest   - The index of the last step that was observed
    %    earliest - The index of the earliest step that has not been forgoten

    properties (Access = private)
        handle;
    end

    methods (Access = public)
        function kalman = UltimateKalmanNative(options)
            if nargin==0
                options = struct();
            end

            if ~isfield(options,'algorithm')
                warning('Setting options.algorithm to Ultimate');
                options.algorithm = 'Ultimate';
            end
            if ~isfield(options,'smoothOnly')
                options.smoothOnly = false;
            end
            if ~isfield(options,'estimateCovariance')
                options.estimateCovariance = true;
            end

            binaryOptions = 0;
            if strcmp(options.algorithm,'Ultimate');     binaryOptions = bitshift(1,0); end
            if strcmp(options.algorithm,'Conventional'); binaryOptions = bitshift(1,1); end
            if strcmp(options.algorithm,'Oddeven');      binaryOptions = bitshift(1,2); end
            if strcmp(options.algorithm,'Associative');  binaryOptions = bitshift(1,3); end
            if binaryOptions==0
                error('The value %s is invalid for options.algorithm',options.algorithm);
            end
            if options.estimateCovariance == false;      binaryOptions = binaryOptions + bitshift(1,16); end
            %if options.smoothOnly         == true;       binaryOptions = binaryOptions + bitshift(1,17); end

            kalman = kalman@handle();
            kalman.handle = ultimatekalmanmex('create',binaryOptions);
            if kalman.handle < 0
                error('no more handles');
            end
        end

        function delete(kalman)
            ultimatekalmanmex('free',kalman.handle);
            kalman.handle = -1; % invalid; not needed but added for safety
        end

        function i = earliest(kalman)
            %EARLIST   The index of the oldest step that has not been
            %          forgotten
            %
            %   kalman.EARLIEST() returns the index of the earliest step
            %   that has not been forgoten. Step numbers start from zero.
            i = ultimatekalmanmex('earliest',kalman.handle);
        end

        function i = latest(kalman)
            %LATEST   The index of the last step that was observed.
            %
            %   kalman.LATEST() returns the index of the latest step that 
            %   was observed.
            i = ultimatekalmanmex('latest',kalman.handle);
        end

        function evolve(kalman,n_i,H_i,F_i,c_i,K_i)
            %EVOLVE   Evolve the state using the given linear recurrence.
            %   kalman.EVOLVE(n_i,H_i,F_i,c_i,K_i) evolves the state using the recurrence
            %                H_i * u_i = F_i * u_{i-1} + c_i + epsilon
            %   where n_i is the dimension of u_i,
            %   C_i is the covariance matrix of the error term epsilon.
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
               ultimatekalmanmex('evolve',kalman.handle,n_i,[],[],[],[],double('X'));
               return;
            end
            l_i = size(F_i,1);                                             % row dimension
            if size(H_i,1) ~= l_i                                          % this allows the user to pass [] for H
                if l_i == n_i
                    H_i = eye(l_i);
                else
                    H_i = [ eye(l_i) zeros(l_i,n_i - l_i)];
                end
            end

            [R,type] = rep(K_i);
            %type
            ultimatekalmanmex('evolve',kalman.handle,n_i,H_i,F_i,c_i,R,double(type));

        end

        function observe(kalman,G_i,o_i,C_i)
            %OBSERVE   Provide observations of the current state.
            %   kalman.OBSERVE(G_i,o_i,C_i) provide observations that
            %   satisfy the linear equation
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
                ultimatekalmanmex('observe',kalman.handle,[],[],[],double('X'));
                return;
            end
            [R,type] = rep(C_i);
            ultimatekalmanmex('observe',kalman.handle,G_i,o_i,R,double(type));

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
            estimate = ultimatekalmanmex('estimate',kalman.handle,s);
            %W        = ultimatekalmanmex('covariance',kalman.handle,s);
            %cov = CovarianceMatrix(W,'W');
            [W,dtype]  = ultimatekalmanmex('covariance',kalman.handle,s);
            cov = CovarianceMatrix(W,char(dtype));
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
 
            ultimatekalmanmex('forget',kalman.handle,s);
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
            ultimatekalmanmex('rollback',kalman.handle,s);
        end
        
        function smooth(kalman)
            %SMOOTH  Compute smooth estimates of all the stored states.
            %
            %   kalman.SMOOTH() computes the smoothed estimated state of all
            %   the steps that are still in memory.
            % 
            %   This method must be called after the last step has been
            %   observed.
            ultimatekalmanmex('smooth',kalman.handle);
        end

        function u = gather(kalman)
            e = kalman.earliest();
            l = kalman.latest();
            N = 0;
            for i=kalman.earliest():kalman.latest()
                N = N + length(kalman.estimate(i));
                 %[range' u(range)]
            end
            %N = kalman.steps{l}.estimateStart + kalman.steps{l}.dimension - 1;
            u = zeros(N,1);

            j = 1;
            for i=kalman.earliest():kalman.latest()
                estimate = kalman.estimate(i);
                range = j:j + length(estimate) - 1;
                u(range) = estimate;
                %[range' u(range)]
                j = j + length(estimate);
                
            end
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
            t = ultimatekalmanmex('perftest',kalman.handle,H,F,c,K_rep,double(K_type),G,o,C_rep,double(C_type),double(count),double(decimation));
        end
        
    end % methods
end