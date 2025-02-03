classdef KalmanExplicitRepresentation < KalmanBase
    % KalmanExplicitRepresentation   A class that stores the evolution and
    %         observation constraints of a Kalman filter/smoother.
    %         This class is meant to be used mostly for Kalman smoothing,
    %         normally not for filtering.
    % 
    % KalmanExplicitRepresentation inherits from KalmanBase.
    %
    % KalmanExplicitRepresentation fields (properties):
    %    nonlinear - A boolean that remembers whether any of the
    %                constraints is nonlinear
    %  
    % KalmanExplicitRepresentation implements these abstract methods:
    %    evolve   - Evolve the state using a linear matrix equation 
    %    observe  - Provide observations of the current state 
    %    smooth   - Compute smooth estimates of all the stored states
    % 
    % KalmanExplicitRepresentation defines additional methods:
    %    initialGuess       - Sets the initial guess for the state of a given step
    %    gatherInitialGuess - Concatenates initial guesses into a column vector 
    %    scatter            - Sets the estimates of all the states from a
    %                         given vector
    %
    % Copyright 2024 Sivan Toledo & Shahaf Gargir, Tel Aviv University.


    properties (Access = public)
        nonlinear;
    end

    methods (Access = protected)
        function rollbackStepWithIndex(kalman,i)
            kalman.current = [];
            kalman.current.dimension = kalman.steps{ i }.dimension;
            kalman.current.step      = kalman.steps{ i }.step;
            kalman.current.H         = kalman.steps{ i }.H;
            kalman.current.F         = kalman.steps{ i }.F;
            kalman.current.c         = kalman.steps{ i }.c;
            kalman.current.K         = kalman.steps{ i }.K;
            kalman.current.l         = kalman.steps{ i }.l;
            kalman.current.estimateStart  = kalman.steps{ i }.estimateStart;
            kalman.current.evolutionStart  = kalman.steps{ i }.evolutionStart;
        end

        function [rows,cols] = size(kalman)
            l = length(kalman.steps);
            cols = kalman.steps{l}.estimateStart + kalman.steps{l}.dimension - 1;
            rows = kalman.steps{l}.observationStart + kalman.steps{l}.m - 1;
        end
    end

    methods (Access = public)
        function kalman = KalmanExplicitRepresentation(options)
            kalman = kalman@KalmanBase(options);
            kalman.nonlinear = false;
            %fprintf('KalmanAccumulator constructor (%s)\n',label);
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
                kalman.current.l = 0;
                kalman.current.estimateStart  = 1;
                kalman.current.evolutionStart = 1;
                return
            end

            ptr_imo = length(kalman.steps);
            kalman.current.step  = kalman.steps{ptr_imo}.step + 1;

            kalman.current.l              = size(c_i,1);
            kalman.current.estimateStart  = kalman.steps{ptr_imo}.estimateStart    + kalman.steps{ptr_imo}.dimension;
            kalman.current.evolutionStart = kalman.steps{ptr_imo}.observationStart + kalman.steps{ptr_imo}.m;

            l_i = kalman.current.l;                                        % row dimension
            if size(H_i,1) ~= l_i                                          % this allows the user to pass [] for H
                if l_i == n_i
                    H_i = eye(l_i);
                else
                    H_i = [ eye(l_i) zeros(l_i,n_i - l_i)];
                end
            end

            kalman.current.H = H_i;
            kalman.current.F = F_i;
            kalman.current.c = c_i;
            kalman.current.K = K_i;     

            if isa(F_i,'function_handle')
                kalman.nonlinear = true;
            end
        end

        function initialGuess(kalman,u_i,s)
            if nargin<3
                kalman.current.initial_guess = u_i;
            else
                kalman.steps{ kalman.indexOfStep(s) }.initial_guess = u_i;
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
            if nargin>=4 && ~isempty(o_i) % we have observations
                kalman.current.G = G_i;
                kalman.current.o = o_i;
                kalman.current.C = C_i;

                kalman.current.m                = size(o_i,1);

                if isa(G_i,'function_handle')
                    kalman.nonlinear = true;
                end

            else
                kalman.current.m                = 0;
            end % of we have observations

            kalman.current.observationStart = kalman.current.evolutionStart + kalman.current.l;

            kalman.steps{ length(kalman.steps)+1 } = kalman.current;
        end

        function u = gatherInitialGuess(kalman)
            l = length(kalman.steps);
            N = kalman.steps{l}.estimateStart + kalman.steps{l}.dimension - 1;
            u = zeros(N,1);
            
            for i=1:l
                range = kalman.steps{i}.estimateStart:kalman.steps{i}.estimateStart + kalman.steps{i}.dimension - 1;
                u(range) = kalman.steps{i}.initial_guess; % step indexes start at zero
                %[range' u(range)]
            end
        end

        function scatter(kalman,u)
            l = length(kalman.steps);

            for i=1:l
                range = kalman.steps{i}.estimateStart:kalman.steps{i}.estimateStart + kalman.steps{i}.dimension - 1;
                kalman.steps{i}.estimatedState = u(range);
            end
        end        

        function [r, J,rows,cols] = residualAndJacobian(kalman,u)
            %fprintf('residualAndJacobian (%s)\n',kalman.label)

            l = length(kalman.steps);

            N = kalman.steps{l}.observationStart + kalman.steps{l}.m - 1;
            r = zeros(N,1);
            if nargout > 1 % we need to return a Jacobian evaluated at u
                J = KalmanJacobian(kalman.options);
            end

            for i=1:l
                range = kalman.steps{i}.estimateStart:kalman.steps{i}.estimateStart + kalman.steps{i}.dimension - 1;
                u_i = u(range);
                %n_i = kalman.steps{i}.dimension;
                l_i = kalman.steps{i}.l;
                m_i = kalman.steps{i}.m;

                if i>1
                    H_i   = kalman.steps{i}.H;
                    F_i   = kalman.steps{i}.F;
                    c_i   = kalman.steps{i}.c;
                    K_i   = kalman.steps{i}.K;
                    if isa(F_i,'function_handle')
                        [F_u_imo, JF] = F_i(u_imo);
                    else
                        F_u_imo = F_i * u_imo;
                        JF = F_i;
                    end
                    range = kalman.steps{i}.evolutionStart:kalman.steps{i}.evolutionStart + kalman.steps{i}.l - 1;
                    r(range) = K_i.weigh( - F_u_imo  + H_i * u_i - c_i);
                    %[-1 i r(range)']
                    if nargout > 1
                        J.evolve(kalman.steps{i}.dimension,K_i.weigh(H_i),K_i.weigh(JF),zeros(l_i,1),CovarianceMatrix(eye(l_i),'W'));
                    end
                else
                    if nargout > 1
                        J.evolve(kalman.steps{i}.dimension); % first step
                    end
                end



                if m_i > 0 % there are observations in this step
                    G_i   = kalman.steps{i}.G;
                    o_i   = kalman.steps{i}.o;
                    C_i   = kalman.steps{i}.C;
                    if isa(G_i,'function_handle')
                        [G_u_i, JG] = G_i(u_i);
                    else
                        G_u_i = G_i * u_i;
                        JG = G_i;
                    end
                    range = kalman.steps{i}.observationStart:kalman.steps{i}.observationStart + kalman.steps{i}.m - 1;
                    r(range) = C_i.weigh( G_u_i - o_i );
                    %[-2 i r(range)']
                    if nargout > 1
                        J.observe(C_i.weigh(JG),zeros(m_i,1),CovarianceMatrix(eye(m_i),'W'));
                    end
                end

                u_imo = u_i;
            end
            if nargout > 1 % we need to return a Jacobian evaluated at u
                cols = kalman.steps{l}.estimateStart + kalman.steps{l}.dimension - 1;
                rows = kalman.steps{l}.observationStart + kalman.steps{l}.m - 1;
                %adapter = KalmanLevenbergMarquardtAdapter(J,rows,cols);
            end
            %norm_dim_r = [norm(r)/sqrt(length(r)) length(r)]
            %title(normr);
            %drawnow;
        end
        
    end % methods
end