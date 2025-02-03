classdef KalmanAssociativeSmoother < KalmanExplicitRepresentation
    % KalmanAssociativeSmoother   A Kalman/RTS (Rauch, Tung, and Striebel) 
    %                             linear smoother. The implementation is
    %                             based on an article that presents a 
    %                             parallel Kalman filter and smoother.
    % 
    % The article is:
    %   Temporal Parallelization of Bayesian Smoothers
    %   by Simo Särkkä, Ángel F. García-Fernández
    %   IEEE Transactions on Automatic Control, 66(1), pages 299-306, 2021
    %   doi 10.1109/TAC.2020.2976316
    %   https://arxiv.org/abs/1905.13002
    % 
    % The implementation was debugged using Python (JAX) code at:
    %   https://github.com/EEA-sensors/sequential-parallelization-examples/tree/main/python/temporal-parallelization-bayes-smoothers
    % 
    % Copyright 2024 Sivan Toledo, Tel Aviv University.

    methods (Access = private)
        function e = buildFilteringElement(kalman,i)
            % we denote by Z the matrix denoted by C in the article,
            % because we already use C for the covariance of the
            % observations.
            n_i  = kalman.steps{i}.dimension;
            step = kalman.steps{i}.step;

            e = struct();

            if step == 0
                return;
            end

            if step == 1

                G_i = kalman.steps{1}.G;
                o_i = kalman.steps{1}.o;
                C_i = kalman.steps{1}.C;
                                
                W_i_G_i = C_i.weigh(G_i);
                W_i_o_i = C_i.weigh(o_i);
                [Q,R] = qr(W_i_G_i,0); 
                m0 = R \ ( Q' * W_i_o_i );
                P0 = inv(R'*R); 

                kalman.steps{1}.estimatedState      = m0;
                kalman.steps{1}.estimatedCovariance = CovarianceMatrix( P0, 'C');
            end

            F_i = kalman.steps{i}.F;
            c_i = kalman.steps{i}.c;
            K_i = kalman.steps{i}.K.explicit;

            if step == 1
                K_i = K_i + F_i * P0 * F_i';
            end

            if ~isfield(kalman.steps{i},'o') % no observations, pass []
                e.Z = K_i;
                if step == 1
                    e.A   = zeros(n_i,n_i);
                    e.b   = m0 + c_i;
               else
                    e.A   = F_i;
                    e.b   = c_i;
                end
                e.e = [];
                e.J = [];
            else % there are observations
                G_i = kalman.steps{i}.G;
                o_i = kalman.steps{i}.o;
                C_i = kalman.steps{i}.C;

                S = G_i * K_i * G_i' + C_i.explicit;
                %G_i_tras_inv_S = G_i' * inv(S)
                G_i_trans_inv_S = (G_i') / S;
                K = K_i * G_i_trans_inv_S ;

                if step == 1
                    e.A   = zeros(n_i,n_i);
                    m1                 = F_i * m0 + c_i;
                    e.b   = m1 + K * ( o_i - G_i * m1 ); 
                    e.Z   = K_i - K * S * K';
               else
                    e.A   = F_i - K * G_i * F_i;
                    e.b   = c_i + K * ( o_i - G_i * c_i );
                    e.Z   = K_i - K * G_i * K_i;
                end
                e.e = F_i' * G_i_trans_inv_S  * ( o_i - G_i * c_i );
                e.J   = F_i' * G_i_trans_inv_S  * G_i * F_i;

            end % of we have observations
        end

        function e = buildSmoothingElement(kalman,i)
            e = struct();

            if i==length(kalman.steps)
                ni = kalman.steps{i}.dimension;
                e.E = zeros(ni,ni);
                e.g = kalman.steps{i}.estimatedState;
                e.L = kalman.steps{i}.estimatedCovariance.explicit();
            else
                x = kalman.steps{i}.estimatedState;
                P = kalman.steps{i}.estimatedCovariance.explicit();
                F = kalman.steps{i+1}.F;
                Q = kalman.steps{i+1}.K.explicit();
                c = kalman.steps{i+1}.c;

                %kalman.steps{i}.E = P * F' * inv(F*P*F'+Q);
                e.E = P * F' / (F*P*F'+Q);
                e.g = x - e.E * (F * x + c);
                e.L = P - e.E * F * P;
            end
        end
    end

    methods (Access = public)
        function kalman = KalmanAssociativeSmoother(options)
            if nargin==0
                options = struct();
            end
            kalman = kalman@KalmanExplicitRepresentation(options);
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

            elements = cell(l);

            for j=1:l
                elements{j} = kalman.buildFilteringElement(j);
            end

            filtered = cummulativeSums(elements, @filteringAssociativeOperation, 2, l, 1);

            i = 1;
            for j=2:l % step{1}.estimatedState and covariance were set in buildFilteringElements(1) above 
                kalman.steps{j}.estimatedState      = filtered{i}.b;
                kalman.steps{j}.estimatedCovariance = CovarianceMatrix( filtered{i}.Z, 'C' );
                i = i+1;
            end

            elements = cell(l);

            for j=1:l
                elements{j} = kalman.buildSmoothingElement(j);
            end

            smoothed = cummulativeSums(elements, @smoothingAssociativeOperation, l, 1, -1);

            i = 2; % TODO!!! do we need to start from 2 (no need to smooth step l)
            for j=l-1:-1:1
                kalman.steps{j}.estimatedState      = smoothed{i}.g;
                kalman.steps{j}.estimatedCovariance = CovarianceMatrix( smoothed{i}.L, 'C' );
                i = i+1;
            end

        end

    end % methods

end

function sums = cummulativeSums(a, f, s, e, stride)
% computes the (inclusive) prefix sums of a cell array 'a'
% from index s to index e with a given stride,
% where the associative operation is f.
%
% The function returns a cell array with the prefix sums. This cell
% array is compact and in order of the sums. Indexes in it do not 
% correspond directly to indexes in a, unless s=stride=1.

    sums = cell(abs(e-s+1)); % could be negative

    i = 1;
    sum = a{s};
    sums{i} = sum;
    i = i+1;

    for j=s+stride:stride:e
        sum = f(sum, a{j});
        sums{i} = sum;
        i = i+1;
    end
end

function sums = cummulativeSumsParallel(a, f, s, e, stride)
% computes the (inclusive) prefix sums of a cell array 'a'
% from index s to index e with a given stride,
% where the associative operation is f.
%
% The function returns a cell array with the prefix sums. This cell
% array is compact and in order of the sums. Indexes in it do not 
% correspond directly to indexes in a, unless s=stride=1.
%
% This implementation can be easily parallelized and has span log_2(n)
% but is not work efficient; it performs \Theta(n log n) work. 
% It does work for any array size.

    sums = cell(1,e-s+1);

    n = (e-s+1)/stride;

    if n <= 10

        i = 1;
        sum = a{s};
        sums{i} = sum;
        i = i+1;

        for j=s+stride:stride:e
            sum = f(sum, a{j});
            sums{i} = sum;
            i = i+1;
        end
    else
        nhalf = floor(n/2);
        m = s + stride*nhalf;

        s1 = cummulativeSums(a, f, s, m, stride);
        s2 = cummulativeSums(a, f, m+stride, e, stride);

        for i=1:length(s2) % in parallel
            s2{i} = f(s1{1,length(s1)},s2{1,i});
        end

        sums = [ s1  s2 ];
    end
end


function sij = filteringAssociativeOperation( si, sj )
    sij = struct();

    ni = length(si.b);

    X = sj.A  / (eye(ni) + si.Z*sj.J);
    Y = si.A' / (eye(ni) + sj.J*si.Z);

    sij.A = X * si.A;
    sij.b = X * (si.b + si.Z*sj.e) + sj.b;
    sij.Z = X * si.Z * sj.A' + sj.Z;
    sij.e = Y * (sj.e - sj.J*si.b) + si.e;
    sij.J = Y * sj.J*si.A + si.J;

end

function sij = smoothingAssociativeOperation( si, sj ) % REVERSED!
    sij = struct();

    sij.E = sj.E * si.E;
    sij.g = sj.E * si.g + sj.g;
    sij.L = sj.E * si.L * sj.E' + sj.L;
end
