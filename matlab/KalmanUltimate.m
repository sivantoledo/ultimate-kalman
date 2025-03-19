classdef KalmanUltimate < KalmanBase
    % KalmanUltimate   An implementation of the Paige-Saunders Kalman
    %                  filter and smoother by Sivan Toledo, Tel Aviv University.
    %
    % KalmanUltimate implements these abstract methods:
    %    evolve   - Evolve the state using a linear matrix equation 
    %    observe  - Provide observations of the current state 
    %    smooth   - Compute smooth estimates of all the stored states
    %
    % Meaningful fields in options structure:
    %    covarianceEstimates - string that speficies the algorithm to
    %                 compute covariance matrices. Valid values are:
    %                 'PaigeSaunders' (default) using orthogonal
    %                     transformations
    %                 'SelInv' using a selective inversion algorithm
    %                 matrices of state estimates (slow if set to true) 
    %
    % For a full documentation of the API, see KalmanBase
    %
    % For a thorough description, see the article:
    %	Sivan Toledo. Algorithm 1051: UltimateKalman, flexible kalman filtering and 
    %   smoothing using orthogonal transformations. ACM Transactions on 
    %   Mathematical Software, 50(4):1-19, 2024.
    %   https://doi.org/10.1145/3699958
    % 
    % Copyright 2020-2024 Sivan Toledo, Tel Aviv University.


    properties (Access = protected)
    end

    methods (Access = protected)
        function rollbackStepWithIndex(kalman,i)
            kalman.current = [];
            kalman.current.dimension = kalman.steps{ i }.dimension;
            kalman.current.step      = kalman.steps{ i }.step;
            kalman.current.Rbar      = kalman.steps{ i }.Rbar;
            kalman.current.ybar      = kalman.steps{ i }.ybar;
        end
    end

    methods (Access = public)
        function kalman = KalmanUltimate(options)
            if nargin==0
                options = struct();
            end
            kalman = kalman@KalmanBase(options);
            kalman.steps = {};
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
                    y = W_i_o_i;
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

            if ~isfield(kalman.options,'covarianceEstimates') || strcmp(kalman.options.covarianceEstimates,'PaigeSaunders')
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

            if isfield(kalman.options,'covarianceEstimates') && strcmp(kalman.options.covarianceEstimates,'SelInv')
            % these formulas are derived from the paper "SelInv—An
            % Algorithm for Selected Inversion of a Sparse Symmetric
            % Matrix" by Lin et al, ACM Trans. on Math Sotware, 2011, 
            % https://doi.org/10.1145/1916461.1916464
            for i=l:-1:1
                if i == l
                    C = kalman.steps{l}.estimatedCovariance.explicit();
                    % the covariance matrix has already been constructed
                    % here.
                else
                    n_i = size(kalman.steps{i}.Rdiag,1);

                    Rdiag    = kalman.steps{i}.Rdiag;
                    Rsupdiag = kalman.steps{i}.Rsupdiag;

                    % RdiagInv = inv(Rdiag);
                    % C = RdiagInv*(eye(n_i) + Rsupdiag*C*Rsupdiag')*RdiagInv'; 
                    C = Rdiag \ ( (eye(n_i) + Rsupdiag*C*Rsupdiag') / (Rdiag') ); 

                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix( C, 'C' );
                end
            end
            end
        end
    end % methods
end