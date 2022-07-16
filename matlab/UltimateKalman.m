classdef UltimateKalman < handle
    % UltimateKalman   An implementation of the Paige-Saunders Kalman
    % filter and smoother by Sivan Toledo, Tel Aviv University.
    %
    % UltimateKalman Methods:
    %    evolve   - Evolve the state using a linear matrix equation
    %    observe  - Provide observations of the current state
    %    filtered - Obtain an estimate of the current state
    %    forget   - Forget the oldest steps to save memory meory
    %    smooth   - Compute smooth estimates of all the stored states
    %    smoothed - Obtain the smoothed estimates of historical states


    properties (Access = private)
        steps;   % old states
        current; % the current state; moved to steps at the end of observe

        k;     % old, now sure what this was
    end

    methods (Access = public)
        function kalman = UltimateKalman()
            kalman = kalman@handle();
            kalman.steps = {};
        end

        %function advance(kalman,dim,timestamp)
            %ADVANCE   Create a new time step with a given state dimension.
            %   kalman.ADVANCE(dimension,timestamp) creates a new time step with
            %   a given dimension. The timestamp can be retrieved by the
            %   user later but is not used by the algorithm. 
            % 
            %   This method must be called before evolve and observe.
        %    kalman.current = []; % clear
        %    if length(kalman.steps) == 0 % starting the first step
        %        kalman.current.step = 0;
        %        kalman.current.local = 1; % local index in steps
        %    else
        %        l = length(kalman.steps);
        %        kalman.current.last = l;
        %        kalman.current.local = l + 1;
        %        kalman.current.step  = kalman.steps{l}.step + 1;
        %    end
        %    kalman.current.dimension = dim;
        %end

        function evolve(kalman,ni,Hi,Fi,ci,Ki)
            %EVOLVE   Evolve the state using the given linear recurrence.
            %   kalman.EVOLVE(H,F,be,C) evolves the state using the recurrence
            %                H * u_i = F * u_{i-1} + be + e
            %   where C is the covariance matrix of the error term e.
            %   The covariance matrix C must be an instance of
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
            % 
            %   This method must be called after advance and before observe.

            kalman.current = []; % clear
            kalman.current.dimension = ni;
            if length(kalman.steps) == 0 % starting the first step
                kalman.current.step = 0;
                kalman.current.local = 1; % local index in steps
                return
            end
                
            l = length(kalman.steps);
            kalman.current.last = l;  % TODO do we really need this? probably not
            kalman.current.local = l + 1;
            kalman.current.step  = kalman.steps{l}.step + 1;

            lastDim = kalman.steps{l}.dimension;

            evolutionCount = size(Fi,1); % row dimension
            if size(Hi,1) ~= evolutionCount % this allows the user to pass [] for H
                if evolutionCount == kalman.current.dimension
                    Hi = eye(evolutionCount);
                else
                    Hi = [ eye(evolutionCount) zeros(evolutionCount,kalman.current.dimension - evolutionCount)];
                end
            end

            WF  = - Ki.weigh(Fi);
            Wbe =   Ki.weigh(ci);
            WI  =   Ki.weigh(Hi);

            l = kalman.current.last;

            if isfield(kalman.steps{l},'Rdiag') && ~isempty(kalman.steps{l}.Rdiag)
                lastRdiagRowDim = size(kalman.steps{l}.Rdiag,1);
                Ki = [ kalman.steps{l}.Rdiag ; WF ];
            else
                lastRdiagRowDim = 0;
                Ki = WF;
            end

            [Q,R] = qr(Ki);
            kalman.steps{l}.Rdiag = R(1:min(size(Ki,1),lastDim),:);
            Z = [ zeros( lastRdiagRowDim, ni) ; WI ];
            Y = Q' * Z;
            kalman.steps{l}.Rsupdiag = Y(1:min(size(Y,1),lastDim),:);

            if (size(Y,1) > lastDim) % we have leftover rows that go into the current Rdiag
                kalman.current.Rdiag = Y(lastDim+1:end,:);
            end

            if isfield(kalman.steps{l},'QTb') && ~isempty(kalman.steps{l}.QTb)
                rhs = [ kalman.steps{l}.QTb ; Wbe ];
            else
                rhs = Wbe;
            end

            v = Q' * rhs;

            kalman.steps{l}.QTb = v(1:min(length(v),lastDim),1);
            if (length(v) > lastDim) % we have leftover rows that go into the current Rdiag
                kalman.current.QTb = v(lastDim+1:end,1);
            else
                kalman.current.QTb = [];
            end

            % not sure if we need to test for the existance of Rdiag
            %if isfield(kalman.current,'Rdiag') && ~isempty(kalman.current.Rdiag) && size(kalman.current.Rdiag,1) == size(kalman.current.Rdiag,1) % Rdiag is square, we can predict the next state
            %    kalman.current.predictedEstimate   = (kalman.current.Rdiag) \ (kalman.current.QTb);
            %    kalman.current.predictedCovariance = CovarianceMatrix(kalman.current.Rdiag, 'W');
            %end
        end

        function observe(kalman,G,bo,C)
            %OBSERVE   Provide observations of the current state.
            %   kalman.OBSERVE(G,bo,C) provide observations that satisfy
            %   the linear equation
            %                bo = G * u_i + e
            %   where C is the covariance matrix of the error term e.
            %   The covariance matrix C must be an instance of
            %   CovarianceMatrix. The matrix G is a standard Matlab
            %   matrix, and bo is a standard column vector. 
            % 
            %   kalman.OBSERVE() tells the algorithm that no
            %   observations are availble of the state of this step.
            % 
            %   This method must be called after advance and evolve.
            d = kalman.current.dimension;
            if nargin<4 || length(bo) == 0 % no observations, pass []
                kalman.current.observationCount = 0;

                if isfield(kalman.current,'Rdiag') && ~isempty(kalman.current.Rdiag)
                    [Q,R] = qr( kalman.current.Rdiag );
                    kalman.current.Rdiag = R;
                    kalman.current.QTb   = Q' * kalman.current.QTb ;
                end
            else % no observations
                kalman.current.observationCount = length(bo);

                WG  = C.weigh(G);
                Wbo = C.weigh(bo);

                if isfield(kalman.current,'Rdiag') && ~isempty(kalman.current.Rdiag)
                    [Q,R] = qr( [ kalman.current.Rdiag ; WG ] , 0 ); % thin QR
                    kalman.current.Rdiag = R;
                    kalman.current.QTb   = Q' * [ kalman.current.QTb ; Wbo ];
                else % no Rdiag yet
                    RdiagRowDim = min(d, size(WG,1));
                    [Q,R] = qr(WG);
                    kalman.current.Rdiag = R(1:RdiagRowDim,1:d);
                    QTb                  = Q' * Wbo;
                    kalman.current.QTb   = QTb(1:RdiagRowDim,1);
                end
            end % of we have observations
            kalman.current.estimatedCovariance = CovarianceMatrix( kalman.current.Rdiag, 'W');

            kalman.steps{ kalman.current.local } = kalman.current;
        end

        function [estimate,cov] = filtered(kalman)
            %FILTERED   Obtain an estimate of the current state and the 
            %           covariance of the estimate.
            %
            %   kalman.FILTERED() returns the estimate of the state of the
            %   last step that was observed (including with an empty observation).
            %   The method also returns the covariance of the estimate. The
            %   covariance is an instance of CovarianceMatrix.
            % 
            %   This method can be called after observe.

            l = length(kalman.steps);
            estimate = kalman.steps{l}.Rdiag \ kalman.steps{l}.QTb;
            cov = kalman.steps{l}.estimatedCovariance;
        end

        function [estimate,cov] = smoothed(kalman,i)
            %SMOOTHED  Obtain an estimate of the smooth estimate of a historical
            %          state and the covariance of the estimate.
            %
            %   kalman.SMOOTHED(i) returns the estimate of state i, which must still be
            %   in memory. It also returns the covariance of the estimate. 
            %   The estimate is only available if the smooth method was called
            %   when step i was in memory.
            % 
            %   This method can be called after smooth.

            first = kalman.steps{1}.step;
            index = i - first + 1;
            estimate = kalman.steps{index}.estimatedState;
            cov = kalman.steps{index}.estimatedCovariance;
        end

        function forget(kalman,howmany)
            %FORGET  Forget the oldest steps to save memory meory.
            %
            %   kalman.FORGET() forgets all but the last step.
            %   
            %   kalman.FORGET(howmany) forgets the requested number of steps
            %   from memory. The steps that are forgoten are the oldest ones.
            %   The method always leaves at least the most recent step in memory. 
            
            l = length(kalman.steps);
            if nargin<2 || howmany >= l
                howmany = l-1;
            end
            kalman.steps = { kalman.steps{howmany+1:l} }; % leave a single step
        end

        function drop(kalman)
            %DROP  Drop all but the last step from meory.
            %
            %   kalman.DROP() drops all but the last step from memory.
            
            l = length(kalman.steps);
            kalman.steps = { kalman.steps{l} }; % leave a single step
        end

        function smooth(kalman)
            %SMOOTH  Compute smooth estimates of all the stored states.
            %
            %   kalman.SMOOTH() computes the smoothed estimated state of all
            %   the steps that are still in memory.

            l = length(kalman.steps);

            v = [];
            for i=l:-1:1
                if length(v)==0 || ~isfield( kalman.steps{i}, 'Rsupdiag' )
                    v = kalman.steps{i}.QTb;
                else
                    v = kalman.steps{i}.QTb - (kalman.steps{i}.Rsupdiag) * v;
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
                    n = size(R,1);
                    m = size(kalman.steps{i}.Rdiag,1);

                    [Q,~] = qr( [ kalman.steps{i}.Rsupdiag ; R ]);

                    QTC = Q' * [ kalman.steps{i}.Rdiag ; zeros(n,size(kalman.steps{i}.Rdiag,2)) ];
                    R = QTC( n+1:n+m , 1:m );

                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix( R, 'W' );
                end
            end
        end

        %============== OLD STUFF BELOW ==============

        function u = old_smooth(obj)
            %obj.k
            %obj.steps{obj.k}
            d = length(obj.steps{obj.k}.QTb);
            for i=obj.k:-1:1
                %i
                %obj.steps{i}
                u(i*d-d+1:i*d,1) = obj.steps{i}.Rdiag \ obj.steps{i}.QTb;
                %[i length(u)]
                if i>1
                    j=i-1;
                    obj.steps{j}.QTb = obj.steps{j}.QTb - obj.steps{j}.Rsupdiag * u(i*d-d+1:i*d,1);
                end
            end
        end

        function invSquareRoots = covInvSquareRoots(obj)
            %obj.k
            %obj.steps{obj.k}
            d = length(obj.steps{obj.k}.QTb);
            Rdiag    = obj.steps{ obj.k }.Rdiag;
            invSquareRoots{ obj.k } = Rdiag;
            for i=obj.k:-1:2
                Rsupdiag_i = obj.steps{ i-1 }.Rsupdiag;
                Rdiag_i    = obj.steps{ i-1 }.Rdiag;
                [ iblock, ~ ] = qr( sparse([ Rsupdiag_i ; Rdiag ]), [ Rdiag_i ; zeros(d,d) ]);

                invSquareRoots{ i-1 } = iblock(d+1:d+d, 1:d);
            end
        end

        function [covs,stddevs] = covariances(obj)
            stddevs = [];
            invSquareRoots = covInvSquareRoots(obj);
            for i=1:obj.k
                isqrt = invSquareRoots{i};

                covsqrt = inv(isqrt);
                cov = covsqrt*covsqrt';
                %full(cov)
                covs{ i } = cov;
                stddevs = [ stddevs ; sqrt(diag(cov)) ];
            end
        end

        function [A,b] = rawLS(obj)
            if isfield(obj.steps{1},'WG')
                A = obj.steps{1}.WG;
                b = obj.steps{1}.Wbo;
            else
                A = [];
                b = [];
            end
            for i=2:obj.k
                [m,n] = size(A);
                d = size(obj.steps{i}.WF,1);
                A = [ A                            zeros(m,d)
                    zeros(d,n-d) -obj.steps{i}.WF obj.steps{i}.WI ];
                b = [ b
                    obj.steps{i}.Wbe ];
                if isfield(obj.steps{i},'WG')
                    l = size(obj.steps{i}.WG,1);
                    A = [ A
                        zeros(l,n) obj.steps{i}.WG ];
                    b = [ b
                        obj.steps{i}.Wbo ];
                end
            end
        end
        function [R,QTb] = triangularLS(obj)
            R   = obj.steps{1}.Rdiag;
            QTb = obj.steps{1}.QTb;
            %size(R)
            for i=2:obj.k-1
                %i
                [m,n] = size(R);
                d = size(obj.steps{i}.Rdiag,1);
                %[m n d]
                R = [ R         [ zeros(m-d,d) ; obj.steps{i-1}.Rsupdiag ]
                    zeros(d,n) obj.steps{i}.Rdiag ];
                QTb = [ QTb
                    obj.steps{i}.QTb ];
            end
            i = obj.k;
            [m,n] = size(R);
            %disp('zzz');
            %i
            %size(obj.steps{i-1}.Rsupdiag,1)
            l=size(obj.steps{i}.Rdiag,1);
            %[m n d]
            %size(R)
            R = [ R         [ zeros(m-d,d) ; obj.steps{i-1}.Rsupdiag ]
                zeros(l,n) obj.steps{i}.Rdiag ];
            QTb = [ QTb
                obj.steps{i}.QTb ];
        end
    end
end