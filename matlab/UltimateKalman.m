classdef UltimateKalman < handle

    properties
        steps;
        first; % global index of the first element in steps
        current;

        k;     % old, now sure what this was
    end

    methods (Access=public)
        function i = getFirst(kalman)
            if length( kalman.steps ) >= 1
                i = kalman.steps{1};
            else
                i = -1;
            end
        end
        function local = getLast(kalman)
            if length( kalman.steps ) >= 1
                local = kalman.steps{ length(kalman.steps) };
            else
                local = -1;
            end
        end

        function local = step(kalman,i)
            % returns the local index of global index i
            local = i - getFirst(kalman);
        end

    end

    methods
        function kalman = UltimateKalman()
            kalman = kalman@handle();
            kalman.steps = {};
            % create the first entry in steps, which represents the initial
            % state
            %kalman.k = 1;
            %kalman.steps{1}.WF = [];
            %kalman.steps{1}.WI = [];
            %kalman.steps{1}.Rdiag = [];
            %kalman.steps{1}.QTb = [];
        end

        function advance(kalman,dim,timestamp)
            %'advance'
            kalman.current = []; % clear
            if length(kalman.steps) == 0 % starting the first step
                kalman.current.step = 0;
                kalman.current.local = 1; % local index in steps
            else
                l = length(kalman.steps);
                kalman.current.last = l;
                kalman.current.local = l + 1;
                kalman.current.step  = kalman.steps{l}.step + 1;
            end


            kalman.current.dimension = dim;
            %kalman.current.timestamp = timestamp;
            %kalman.current
        end

        function evolve(kalman,H,F,be,C)
            %'evolve'
            d = kalman.current.dimension;
            evolutionCount = size(F,1); % row dimension
            if size(H,1) ~= evolutionCount % this allows the user to pass [] for H
                if evolutionCount == kalman.current.dimension
                    H = eye(evolutionCount);
                else
                    H = [ eye(evolutionCount) zeros(evolutionCount,kalman.current.dimension - evolutionCount)];
                end
            end

            WF  = - C.weigh(F);
            Wbe =   C.weigh(be);
            WI  =   C.weigh(H);

            l = kalman.current.last;

            if isfield(kalman.steps{l},'Rdiag') && ~isempty(kalman.steps{l}.Rdiag)
                lastRdiagRowDim = size(kalman.steps{l}.Rdiag,1);
                K = [ kalman.steps{l}.Rdiag ; WF ];
            else
                lastRdiagRowDim = 0;
                K = WF;
            end

            [Q,R] = qr(K);
            kalman.steps{l}.Rdiag = R(1:min(size(K,1),kalman.steps{l}.dimension),:);
            Z = [ zeros( lastRdiagRowDim, d) ; WI ];
            Y = Q' * Z;
            kalman.steps{l}.Rsupdiag = Y(1:min(size(Y,1),kalman.steps{l}.dimension),:);

            if (size(Y,1) > kalman.steps{l}.dimension) % we have leftover rows that go into the current Rdiag
                kalman.current.Rdiag = Y(kalman.steps{l}.dimension+1:end,:);
            end

            if isfield(kalman.steps{l},'QTb') && ~isempty(kalman.steps{l}.QTb)
                rhs = [ kalman.steps{l}.QTb ; Wbe ];
            else
                rhs = Wbe;
            end

            v = Q' * rhs;

            kalman.steps{l}.QTb = v(1:min(length(v),kalman.steps{l}.dimension),1);
            if (length(v) > kalman.steps{l}.dimension) % we have leftover rows that go into the current Rdiag
                kalman.current.QTb = v(kalman.steps{l}.dimension+1:end,1);
            else
                kalman.current.QTb = [];
            end
        end

        function observe(kalman,G,bo,C)
            %'observe'
            d = kalman.current.dimension;
            if length(bo) == 0 % no observations, pass []
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
            l = length(kalman.steps);
            %step = kalman.steps{ kalman.current.local };
            estimate = kalman.steps{l}.Rdiag \ kalman.steps{l}.QTb;
            cov = kalman.steps{l}.estimatedCovariance;
        end

        function [estimate,cov] = smoothed(kalman,i)
            first = kalman.steps{1}.step;
            index = i - first + 1;
            estimate = kalman.steps{index}.estimatedState;
            cov = kalman.steps{index}.estimatedCovariance;
        end

        function drop(kalman)
            l = length(kalman.steps);
            kalman.steps = { kalman.steps{l} }; % leave a single step
        end

        function smooth(kalman)
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