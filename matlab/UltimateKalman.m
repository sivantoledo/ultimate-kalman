classdef UltimateKalman < handle
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
    %    forget   - Forget the oldest steps to save memory memory 
    %    rollback - Roll back the object to an earlier step 
    %    latest   - The index of the last step that was observed 
    %    earliest - The index of the earliest step that has not been forgoten

    properties (Access = protected)
        steps;   % old states
        current; % the current state; moved to steps at the end of observe
    end

    methods (Access = public)
        function kalman = UltimateKalman()
            kalman = kalman@handle();
            kalman.steps = {};
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

        function evolve(kalman,n_i,H_i,F_i,c_i,K_i)
            %EVOLVE   Evolve the state using the given linear recurrence.
            %   kalman.EVOLVE(n_i,H_i,F_i,c_i,K_i) evolves the state using the recurrence
            %                H_i * u_i = F_i * u_{i-1} + c_i + epsilon
            %   where n_i is the dimension of u_i,
            %   C_i is the covariance matrix of the error term epsilon.
            %   The covariance matrix C_i must be an instance of
            %   CovarianceMatrix. The matrices H_i and F_i are standard Matlab
            %   matrices, and be is a standard column vector.
            % 
            %   In the first step, this method does nothing; the first state is
            %   not an evolution of an earlier state. You can omit
            %   arguments or provide empty arguments in the first step.
            %
            %   The argument H can be empty, []. If it is, H is taken to be
            %   an identity matrix (possibly rectangular, if the dimension of
            %   this state is larger than the dimension of the previous state).

            kalman.current = []; % clear
            kalman.current.dimension = n_i;
            if isempty(kalman.steps) % starting the first step
                kalman.current.step = 0;
                %kalman.current.local = 1; % local index in steps
                return
            end
                
            ptr_imo = length(kalman.steps);                                % pointer to row i-1 in the cell array
            %kalman.current.last = ptr_imo;  % TODO do we really need this? probably not
            %kalman.current.local = ptr_imo + 1;
            kalman.current.step  = kalman.steps{ptr_imo}.step + 1;

            n_imo = kalman.steps{ptr_imo}.dimension;                             % n_{i-1}

            l_i = size(F_i,1);                                             % row dimension
            if size(H_i,1) ~= l_i                                          % this allows the user to pass [] for H
                if l_i == n_i
                    H_i = eye(l_i);
                else
                    H_i = [ eye(l_i) zeros(l_i,n_i - l_i)];
                end
            end

            V_i_F_i = - K_i.weigh(F_i);
            V_i_c_i =   K_i.weigh(c_i);
            V_i_H_i =   K_i.weigh(H_i);

            %ptr_imo = kalman.current.last;

            % we denote by z_i the row dimension of Rtilde_{i-1,i-1}
            if isfield(kalman.steps{ptr_imo},'Rdiag') ...
               && ~isempty(kalman.steps{ptr_imo}.Rdiag)                    % coverage tested in ... ?
                z_i = size(kalman.steps{ptr_imo}.Rdiag,1);    
                A = [ kalman.steps{ptr_imo}.Rdiag ; V_i_F_i ];
                B = [ zeros( z_i, n_i)            ; V_i_H_i ];
                y = [ kalman.steps{ptr_imo}.y ; V_i_c_i ];
            else
                %z_i = 0;
                A = V_i_F_i;
                B = V_i_H_i;
                y = V_i_c_i;
            end

            [Q,R] = qr(A);
            B = Q' * B;
            y = Q' * y;
            kalman.steps{ptr_imo}.Rdiag    = R(1:min(size(A,1),n_imo),:);
            kalman.steps{ptr_imo}.Rsupdiag = B(1:min(size(B,1),n_imo),:);
            kalman.steps{ptr_imo}.y        = y(1:min(length(y),n_imo),1);

            % block row i-1 is now sealed

            if (size(B,1) > n_imo)                                         % we have leftover rows that go into the Rbar
                kalman.current.Rbar = B(n_imo+1:end,:);
                kalman.current.ybar = y(n_imo+1:end,1);
            end

            %if kalman.current.step > 395 && kalman.current.step < 405
            %fprintf("ev %d\n",kalman.current.step);
            %size(A)
            %size(B)
            %end
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
            n_i = kalman.current.dimension;
            if nargin<4 || isempty(o_i) % no observations, pass []
                %m_i = 0;
                if isfield(kalman.current,'Rbar') && ~isempty(kalman.current.Rbar)
                    A = kalman.current.Rbar;
                    y = kalman.current.ybar;

                    %[Q,R] = qr( kalman.current.Rbar );
                    %kalman.current.Rdiag = R;
                    %kalman.current.y     = Q' * kalman.current.ybar ;
                else
                    A = [];
                    y = [];
                end
            else % no observations
                %m_i = length(o_i);

                W_i_G_i = C_i.weigh(G_i);
                W_i_o_i = C_i.weigh(o_i);

                if isfield(kalman.current,'Rbar') && ~isempty(kalman.current.Rbar)
                    A = [ kalman.current.Rbar ; W_i_G_i ];
                    y = [ kalman.current.ybar ; W_i_o_i ];

                    %[Q,R] = qr( [ kalman.current.Rbar ; W_i_G_i ] , 0 ); % thin QR
                    %kalman.current.Rdiag = R;
                    %kalman.current.y   = Q' * [ kalman.current.ybar ; W_i_o_i ];
                else % no Rbar yet
                    A = W_i_G_i;
                    y =  W_i_o_i;
                    %RdiagRowDim = min(n_i, size(W_i_G_i,1));
                    %[Q,R] = qr(W_i_G_i);
                    %kalman.current.Rdiag = R(1:RdiagRowDim,1:n_i);
                    %y                  = Q' * W_i_o_i;
                    %kalman.current.y   = y(1:RdiagRowDim,1);
                end
            end % of we have observations

            %if kalman.current.step > 395 && kalman.current.step < 405
            %fprintf("ob %d\n",kalman.current.step);
            %size(A)
            %end

            %'observe'
            %A
            %y

            if ~isempty(A) 
                if size(A,1) >= size(A,2)
                  [Q,R] = qr(A,0);
                  y = Q'*y;
                  n = min(n_i,size(A,1));
                  %n_i
                  %R
                  %y
                  kalman.current.Rdiag = R(1:n,:);
                  kalman.current.y   = y(1:n,1);
                else
                  kalman.current.Rdiag = A;
                  kalman.current.y   = y;
                end
            end
            

            if isfield(kalman.current,'Rdiag') && size(kalman.current.Rdiag,1) == kalman.current.dimension
                kalman.current.estimatedState      = kalman.current.Rdiag \ kalman.current.y;
                %Rdiag = kalman.current.Rdiag
                kalman.current.estimatedCovariance = CovarianceMatrix( kalman.current.Rdiag, 'W');
            end

            kalman.steps{ length(kalman.steps)+1 } = kalman.current;
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
            l = length(kalman.steps);
            latest   = kalman.steps{l}.step;
            earliest = kalman.steps{1}.step;
            if nargin < 2
                s = latest;
            end
            if s < earliest || s > latest 
                warning('cannot provide an estimate, too old or in the future')
                return
            end
            ptr_s = s - earliest + 1;

            if isfield(kalman.steps{ptr_s},'estimatedState') && ~isempty(kalman.steps{ptr_s}.estimatedState)
                estimate = kalman.steps{ptr_s}.estimatedState;
                cov      = kalman.steps{ptr_s}.estimatedCovariance;
            else
                estimate = NaN * zeros(kalman.steps{ptr_s}.dimension,1);
                cov      = CovarianceMatrix(NaN*eye(kalman.steps{ptr_s}.dimension),'W');
                warning('state is currently underdetermined');
            end
        end

        function forget(kalman,s)
            %FORGET  Forget the oldest steps to save memory meory.
            %   
            %   kalman.FORGET(s) forgets steps s and older. 
            %
            %   kalman.FORGET()  forgets all but the last step.
            
            l = length(kalman.steps);
            earliest = kalman.steps{ 1 }.step;
            if nargin<2
                ptr_s = l-1;
            else
                ptr_s = s - earliest + 1;
            end
            kalman.steps = { kalman.steps{ptr_s+1:l} };
        end

        function rollback(kalman,s)
            %ROLLBACK  Rolls the object back to just after the evolution
            %          to step s (but before the observation in step s).
            %
            %   kalman.ROLLBACK(s) rolls back the filter. Steps later than
            %   s are discarded, and so is the observation of step s, if
            %   there was any.
            % 
            %   This method allows you to roll back the object after it was
            %   used for prediction of future state that have not been
            %   observed yet. Typically steps s+1 and higher do not have
            %   any observation at the time of the rollback.
            
            l = length(kalman.steps);
            earliest = kalman.steps{ 1 }.step;
            ptr_s = s - earliest + 1;
            if ptr_s > l || ptr_s < 1
                warning('cannot roll back to this state (too old or future');
            else
                if true
                kalman.current = [];
                kalman.current.dimension = kalman.steps{ ptr_s }.dimension;
                kalman.current.step      = kalman.steps{ ptr_s }.step;
                kalman.current.Rbar      = kalman.steps{ ptr_s }.Rbar;
                kalman.current.ybar      = kalman.steps{ ptr_s }.ybar;
                %kalman.current = kalman.steps{ ptr_s };
                kalman.steps = { kalman.steps{1:ptr_s-1} };
                else
                kalman.steps = { kalman.steps{1:ptr_s} };
                end
                %kalman.current = rmfield(kalman.current,{'Rdiag','Rsupdiag','y','estimatedState','estimatedCovariance'});
                %ptr_s
                %kalman
                %kalman.steps
                %kalman.current
            end
        end
        
        function smooth(kalman)
            %SMOOTH  Compute smooth estimates of all the stored states.
            %
            %   kalman.SMOOTH() computes the smoothed estimated state of all
            %   the steps that are still in memory.
            % 
            %   This method must be called after the last step has been
            %   observed.

            l = length(kalman.steps);

            v = [];
            for i=l:-1:1
                if i == l
                    v = kalman.steps{i}.y;
                else
                    v = kalman.steps{i}.y - (kalman.steps{i}.Rsupdiag) * v;
                end

                v = (kalman.steps{i}.Rdiag) \ v;

                kalman.steps{i}.estimatedState = v;
            end

            R = [];
            for i=l:-1:1
                if i == l
                    R = kalman.steps{i}.Rdiag;
                    % the covariance matrix has already been constructed
                    % here.
                else
                    n_ipo   = size(R,1);
                    n_i = size(kalman.steps{i}.Rdiag,1);

                    [Q,~] = qr( [ kalman.steps{i}.Rsupdiag ; R ]);

                    S = Q' * [ kalman.steps{i}.Rdiag ; zeros(n_ipo,size(kalman.steps{i}.Rdiag,2)) ];
                    R = S( n_ipo+1:n_ipo+n_i , 1:n_i );

                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix( R, 'W' );
                end
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

    end % methods
end