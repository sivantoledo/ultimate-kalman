classdef KalmanSparse < KalmanBase
    % KalmanSparse   A linear Kalman filter and RTS (Rauch, Tung, and Striebel) 
    %                smoother. The implementation is based on representing
    %                the problem as a sparse least squares problem. The
    %                construction of the sparse matrix is efficient.
    %                Filtering (obtaining an estimate in every step) is
    %                correct but inefficient because it solves a large
    %                least squares problem in every step. Smoothing is
    %                efficient (constructs and factors the matrix once) but
    %                obtaining the covariance matrices is inefficient, and
    %                done through the construction of the global covariance
    %                matrix (not just its diagonal blocks).
    % 
    %                This class can be used to test the correctness and
    %                accuracy of other linear Kalman filters and smoothers.
    %
    % KalmanSparse inherits from KalmanBase.
    %
    % KalmanSparse fields (properties): internal only.
    %  
    % KalmanSparse implements these abstract methods:
    %    evolve   - Evolve the state using a linear matrix equation 
    %    observe  - Provide observations of the current state 
    %    smooth   - Compute smooth estimates of all the stored states
    %
    % Meaningful fields in options structure:
    %    smoothOnly - Boolean, whether to only smooth or also support filtering 
    %                 (filtering is inefficient)
    %    estimateCovariance - Boolean, whether to estimate the covariance
    %                 matrices of state estimates (slow if set to true) 
    % 
    % Copyright 2024 Sivan Toledo and Shahaf Gargir, Tel Aviv University.

    properties (Access = protected)
        %smoothOnly;
        
        i;
        j;
        v;
        ijvLength;

        rhs;
        rhsLength;
    end

    methods (Access = protected)
        function rollbackStepWithIndex(kalman,i)
            kalman.current = [];
            kalman.current.dimension   = kalman.steps{ i }.dimension;
            kalman.current.step        = kalman.steps{ i }.step;
            kalman.current.start       = kalman.steps{ i }.start;
            kalman.current.end         = kalman.steps{ i }.end;
            kalman.current.evolve_start= kalman.steps{ i }.evolve_start;
            kalman.current.evolve_end  = kalman.steps{ i }.evolve_end;
            kalman.current.ijv_e_start = kalman.steps{ i }.ijv_e_start;
            kalman.current.ijv_e_end   = kalman.steps{ i }.ijv_e_end;

            kalman.ijvLength           = kalman.steps{ i }.ijv_e_end;
            kalman.rhsLength           = kalman.steps{ i }.evolve_end;
        end
    end


    methods (Access = public)
        function kalman = KalmanSparse(options)
            if nargin==0
                options = struct();
            end
            kalman = kalman@KalmanBase(options);

            if ~isfield(options,'smoothOnly')
                 kalman.options.smoothOnly = false;
            end
            if ~isfield(options,'estimateCovariance')
                 kalman.options.estimateCovariance = true;
            end

            if kalman.options.smoothOnly == false
                warning('filtering is slow in KalmanSparse, consider setting options.smoothOnly=true')
            end

            if kalman.options.estimateCovariance == true
                warning('estimate covariance matrices is slow in KalmanSparse, consider setting options.estimateCovariance=false')
            end

            kalman.i = zeros(1,1);
            kalman.j = zeros(1,1);
            kalman.v = zeros(1,1);
            kalman.ijvLength = 0;

            kalman.rhs = zeros(1,1);
            kalman.rhsLength = 0;
        end

        % no need for a destructor

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

            kalman.current = []; % clear
            kalman.current.dimension = n_i;
            if isempty(kalman.steps) % starting the first step
                kalman.current.step         = 0;

                kalman.current.start        = 1;   % column index in A, index in u
                kalman.current.end          = n_i; % column index in A
                kalman.current.evolve_start = 1; % row index in A
                kalman.current.evolve_end   = 0; % row index in A
                kalman.current.ijv_e_start  = 1; % row index in ijv representation
                kalman.current.ijv_e_end    = 0; % row index in ijv representation

                %disp('evolve')
                %kalman.current

                return
            end
                
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

            ptr_imo = length(kalman.steps);                                % pointer to row i-1 in the cell array
            kalman.current.step  = kalman.steps{ptr_imo}.step + 1;

            n_imo = kalman.steps{ptr_imo}.dimension;                       % n_{i-1}
            l_i = size(F_i,1);
            
            block = [ V_i_F_i V_i_H_i ];
            [i, j, v] = find(block);

            kalman.current.start        = kalman.steps{ptr_imo}.end         + 1;
            kalman.current.evolve_start = kalman.steps{ptr_imo}.observe_end + 1;
            kalman.current.ijv_e_start  = kalman.steps{ptr_imo}.ijv_o_end + 1;

            kalman.current.end          = kalman.current.start + n_i             - 1;
            kalman.current.evolve_end   = kalman.current.evolve_start + l_i      - 1;
            kalman.current.ijv_e_end    = kalman.current.ijv_e_start + length(v) - 1;

            % shift block indices to global indices
            
            i = i + kalman.current.evolve_start - 1; 
            j = j + kalman.steps{ptr_imo}.start - 1; % the F_i block belongs to the columns of the previous step
            kalman.ijvAppend(i,j,v);

            kalman.rhsAppend(V_i_c_i);

            % compute a predicted estimate

            if ~kalman.options.smoothOnly
                A = sparse( kalman.i(1:kalman.ijvLength), ...
                    kalman.j(1:kalman.ijvLength), ...
                    kalman.v(1:kalman.ijvLength), ...
                    kalman.current.evolve_end,    ...
                    kalman.current.end);
                b = kalman.rhs( 1:kalman.rhsLength );

                [Qtb,R] = qr(A,b);
                u = R \ Qtb;

                kalman.current.estimatedState      = u(kalman.current.start:kalman.current.end);
                kalman.current.estimatedCovariance = CovarianceMatrix( full(R(kalman.current.start:kalman.current.end,kalman.current.start:kalman.current.end)), 'W');
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
            kalman.current.observe_start = kalman.current.evolve_end + 1;
            kalman.current.ijv_o_start   = kalman.current.ijv_e_end  + 1; 

            if nargin<4 || isempty(o_i) % no observations, pass []
                kalman.current.observe_end   = kalman.current.evolve_end; 
                kalman.current.ijv_o_end     = kalman.current.ijv_e_end ; 
            else % no observations
                W_i_G_i = C_i.weigh(G_i);
                W_i_o_i = C_i.weigh(o_i);

                m_i = size(G_i,1);

                [i, j, v] = find(W_i_G_i);
                % shift block indices to global indices
                i = i + kalman.current.observe_start - 1; 
                j = j + kalman.current.start         - 1;

                kalman.current.observe_end   = kalman.current.evolve_end + m_i; 
                kalman.current.ijv_o_end     = kalman.current.ijv_e_end + length(v); 

                kalman.ijvAppend(i,j,v);

                kalman.rhsAppend(W_i_o_i);

            end % of we have observations

            if ~kalman.options.smoothOnly
                A = sparse( kalman.i(1:kalman.ijvLength), ...
                    kalman.j(1:kalman.ijvLength), ...
                    kalman.v(1:kalman.ijvLength), ...
                    kalman.current.observe_end,   ...
                    kalman.current.end);
                b = kalman.rhs( 1:kalman.rhsLength );

                [Qtb,R] = qr(A,b);
                u = R \ Qtb;

                kalman.current.estimatedState      = u(kalman.current.start:kalman.current.end);
                kalman.current.estimatedCovariance = CovarianceMatrix( full(R(kalman.current.start:kalman.current.end,kalman.current.start:kalman.current.end)), 'W');
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

            %warning('unimplemented yet')
            %return;

            l = length(kalman.steps);

            A = sparse( kalman.i(1:kalman.ijvLength), ...
                        kalman.j(1:kalman.ijvLength), ...
                        kalman.v(1:kalman.ijvLength), ...
                        kalman.current.observe_end,                    ...
                        kalman.current.end); 
            b = kalman.rhs( 1:kalman.rhsLength );

            [Qtb,R] = qr(A,b);
            u = R \ Qtb;

            if kalman.options.estimateCovariance

                % the following step is slow for large matrices; we invert a dense
                % matrix to get the entire covariance matrix, not just the
                % diagonal blocks
                cov_u = inv(R' * R);

                for i=l:-1:1
                    kalman.steps{i}.estimatedState      = u(kalman.steps{i}.start:kalman.steps{i}.end);
                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix( full(cov_u(kalman.steps{i}.start:kalman.steps{i}.end,kalman.steps{i}.start:kalman.steps{i}.end)), 'C');
                end
            else
                for i=l:-1:1
                    kalman.steps{i}.estimatedState      = u(kalman.steps{i}.start:kalman.steps{i}.end);
                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix( zeros(kalman.steps{i}.dimension,kalman.steps{i}.dimension), 'C');
                end
            end

            return;
        end

        function A = getA(kalman)
            l = length(kalman.steps);

            A = sparse( kalman.i(1:kalman.ijvLength), ...
                        kalman.j(1:kalman.ijvLength), ...
                        kalman.v(1:kalman.ijvLength), ...
                        kalman.current.observe_end,                    ...
                        kalman.current.end); 

        end

    end % methods

    methods (Access = private)
        function ijvAppend(kalman,ii,jj,vv)
            r = length(vv);
            if kalman.ijvLength + r > length(kalman.v)
                kalman.ijvExtend();
            end
            kalman.i( kalman.ijvLength+1:kalman.ijvLength+r, 1 ) = ii;
            kalman.j( kalman.ijvLength+1:kalman.ijvLength+r, 1 ) = jj;
            kalman.v( kalman.ijvLength+1:kalman.ijvLength+r, 1 ) = vv;

            kalman.ijvLength = kalman.ijvLength + r;
        end

        function ijvExtend(kalman)
            ii = kalman.i;
            jj = kalman.j;
            vv = kalman.v;
            
            r = length(vv);

            kalman.i = zeros(2*r,1);
            kalman.j = zeros(2*r,1);
            kalman.v = zeros(2*r,1);

            kalman.i(1:r,1) = ii;
            kalman.j(1:r,1) = jj;
            kalman.v(1:r,1) = vv;
        end

        function rhsAppend(kalman,rr)
            r = length(rr);
            if kalman.rhsLength + r > length(kalman.rhs)
                kalman.rhsExtend();
            end
            kalman.rhs( kalman.rhsLength+1:kalman.rhsLength+r, 1 ) = rr;

            kalman.rhsLength = kalman.rhsLength + r;
        end

        function rhsExtend(kalman)
            rr = kalman.rhs;
            
            r = length(rr);

            kalman.rhs = zeros(2*r,1);

            kalman.rhs(1:r,1) = rr;
        end

    end
end