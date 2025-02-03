classdef (Abstract) KalmanBase < handle
    % KalmanBase  Base class for Kalman filters and smoothers implemented 
    %             in Matlab (there are separate interface classes for the
    %             native and Java implementations).
    %
    % KalmanBase fields (properties):
    %    steps    - A cell array of structures that represent the time steps
    %    current  - Represents the current time step
    %    options  - Structure containing options
    %  
    % KalmanBase abstract methods:
    %    evolve   - Evolve the state using a linear matrix equation 
    %    observe  - Provide observations of the current state 
    %    smooth   - Compute smooth estimates of all the stored states
    % 
    % KalmanBase methods
    %    estimate - Return the most up to date estimate of a state vector
    %               and its covariance matrix
    %    forget   - Forget the oldest steps to save memory meory 
    %    rollback - Roll back the object to an earlier step 
    %    latest   - The index of the last step that was observed 
    %    earliest - The index of the earliest step that has not been forgoten
    %    gather   - return all the estimates concatenated into a column
    %               vector
    %    pertest  - measure the performance of Kalman filtering or
    %               smoothing
    %
    % Subclasses must define a rollbackStepWithIndex method that should reconstruct
    % the current structure to the point just after the call to evolve() on
    % a given step in the sequence.
    %
    % Copyright 2024 Sivan Toledo, Tel Aviv University.

    properties (Access = protected)
        steps;   % must have at least a member 'step'
        current; % the current state; moved to steps at the end of observe
        options;
    end

    methods (Abstract, Access = public)
        evolve(kalman,n_i,H_i,F_i,c_i,K_i)
        observe(kalman,G_i,o_i,C_i)
        smooth(kalman)
    end

    methods (Abstract, Access = protected)
        rollbackStepWithIndex(kalman,i)
    end

    methods (Access = protected)
        function i = indexOfStep(kalman,s)
            if isempty(kalman.steps)
                error('invalid internal state');
            end

            l = length(kalman.steps);
            earliest = kalman.steps{ 1 }.step;
            i = s - earliest + 1;
            if i > l || i < 1
                error('invalid internal state');
            end
        end
    end


    methods (Access = public)
        function kalman = KalmanBase(options)
            kalman = kalman@handle();
            kalman.steps = {};
            kalman.options = options;
        end

        % no need for a destructor

        function i = earliest(kalman)
            %EARLIST   The index of the oldest step that has not been
            %          forgotten
            %
            %   kalman.EARLIEST() returns the index of the earliest step
            %   that has not been forgoten. Step numbers start from zero.
            if isempty(kalman.steps)
                i = -1;
            else
                i = kalman.steps{ 1 }.step;
            end
        end

        function i = latest(kalman)
            %LATEST   The index of the last step that was observed.
            %
            %   kalman.LATEST() returns the index of the latest step that 
            %   was observed.
            if isempty(kalman.steps)
                i = -1;
            else
                i = kalman.steps{ length(kalman.steps) }.step;
            end
        end

        function forget(kalman,s)
            %FORGET  Forget the oldest steps to save memory meory.
            %   
            %   kalman.FORGET(s) forgets steps s and older. 
            %
            %   kalman.FORGET()  forgets all but the last step.
            
            %l = length(kalman.steps);
            % earliest = kalman.steps{ 1 }.step;
            % if nargin<2
            %     ptr_s = l-1;
            % else
            %     ptr_s = s - earliest + 1;
            % end
            l = length(kalman.steps);
            ptr_s = kalman.indexOfStep(s);
            kalman.steps = { kalman.steps{ptr_s+1:l} };
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
            
            % l = length(kalman.steps);
            % earliest = kalman.steps{ 1 }.step;
            % ptr_s = s - earliest + 1;
            % if ptr_s > l || ptr_s < 1
            %     warning('cannot roll back to this state (too old or future');
            % else
            ptr_s = kalman.indexOfStep(s);
                kalman.rollbackStepWithIndex(ptr_s); % abstract method to return this step to just after evolve
                kalman.steps = { kalman.steps{1:ptr_s-1} }; % this deletes the step we rolled back to!
            % end
        end

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
            
            if nargin<2
                s = kalman.latest();
            end
            
            ptr_s = kalman.indexOfStep(s);

            if isfield(kalman.steps{ptr_s},'estimatedState') && ~isempty(kalman.steps{ptr_s}.estimatedState)
                estimate = kalman.steps{ptr_s}.estimatedState;
                cov      = kalman.steps{ptr_s}.estimatedCovariance;
            else
                estimate = NaN * zeros(kalman.steps{ptr_s}.dimension,1);
                cov      = CovarianceMatrix(NaN*eye(kalman.steps{ptr_s}.dimension),'W');
                warning('state is currently underdetermined');
            end
        end

        function u = gather(kalman)
            l = length(kalman.steps);
            N = 0;
            for i=1:l
                N = N + kalman.steps{i}.dimension;
                 %[range' u(range)]
            end
            %N = kalman.steps{l}.estimateStart + kalman.steps{l}.dimension - 1;
            u = zeros(N,1);

            j = 1;
            for i=1:l
                range = j:j + kalman.steps{i}.dimension - 1;
                u(range) = kalman.estimate(kalman.steps{i}.step);
                %[range' u(range)]
                j = j + kalman.steps{i}.dimension;
                
            end
        end
        
        function t = perftest(kalman,H,F,c,K,G,o,C,count,decimation)
            % Stress testing
            n = size(G,2);
            m = size(G,1);

            %count
            %decimation
            %floor(count/decimation)

            t = zeros(floor(count/decimation),1);

            tstart = tic;
            
            j = 1;

            for i=1:count
                evolve(kalman,n,H,F,c,K);
                observe(kalman,G,o,C);
                [e,cov] = estimate(kalman);
                forget(kalman);
                if rem(i,decimation)==0
                    telapsed = toc(tstart);
                    t(j) = telapsed/decimation;
                    j = j+1;
                    tstart = tic;
                end
            end
        end

        function t = perftestSmoothing(kalman,H,F,c,K,G,o,C,count)
            % Stress testing
            n = size(G,2);

            t = NaN * zeros(2,1);

            tstart = tic;
            
            for i=1:count
                evolve(kalman,n,H,F,c,K);
                observe(kalman,G,o,C);
            end
            t(2,1) = toc(tstart);
            smooth(kalman);

            t(1,1) = toc(tstart);
        end

        
    end % methods
end