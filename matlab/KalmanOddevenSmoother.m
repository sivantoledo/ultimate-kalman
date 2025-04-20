classdef KalmanOddevenSmoother < KalmanExplicitRepresentation
    % KalmanOddevenSmoother   An implementation of the Gargir-Toledo
    %                         parallel Kalman smoother.
    %
    % KalmanOddEvenSmoother implements these abstract methods:
    %    smooth   - Compute smooth estimates of all the stored states
    %
    % KalmanOddEvenSmoother derives from KalmanExplicitRepresentation,
    % which implements the evolve and observe methods.
    %
    % For a full documentation of the API, see KalmanBase.
    %
    % For a thorough description, see the article:
    %	Shahaf Gargir and Sivan Toledo. A Parallel-in-time Kalman smoothing
    %   using orthogonal transformations. In Proceedings of the IEEE
    %   International Parallel and Distributed Processing Symposium (IPDPS), 2025.
    %
    % Copyright 2020-2024 Shahaf Gargir and Sivan Toledo, Tel Aviv University.

    methods (Access = private)
        function parallel_smooth (kalman,indices)
            %SMOOTH  Compute smooth estimates of all the stored states.
            %
            %   kalman.SMOOTH() computes the smoothed estimated state of all
            %   the steps that are still in memory.
            %
            %   This method must be called after the last step has been
            %   observed.

            if(nargin == 1)
                indices = 1:length(kalman.steps);
            end

            %indices

            if (length(indices) == 1)
                singleStep = kalman.steps{indices(1)};
                % singleStep.C should be upper-triangular at this point
                R = singleStep.C;
                singleStep.estimatedState = R \ singleStep.o;
                %fprintf('Xs   set covariance for step %d;\n',indices(1));
                singleStep.S = R \ ( eye(singleStep.dimension) / (R'));
                singleStep.estimatedCovariance = CovarianceMatrix(singleStep.S,'C');
                %fprintf('.   set covariance for step %d\n',indices(1));
                kalman.steps{indices(1)} = singleStep;
                return
            end

            l = length(indices);

            for j=0:2:(l - 1)
                i = indices(j + 1);

                if (j + 1 == l) % Last index
                    C_i = kalman.steps{i}.C;
                    [Q,R_tilde] = qr(C_i);
                    kalman.steps{i}.R_tilde = R_tilde;

                    kalman.steps{i}.o = Q' * kalman.steps{i}.o;
                else
                    i_p_1 = indices(j + 2);
                    C_i = kalman.steps{i}.C;
                    B_i_p_1 = kalman.steps{i_p_1}.B;
                    D_i_p_1 = kalman.steps{i_p_1}.D;

                    mm = size(C_i,1);
                    %ll = size(B_i_p_1,1);
                    nn = size(B_i_p_1,2); % == size(C_i,2)

                    [Q,R_tilde] = qr([C_i ; B_i_p_1]);
                    ll = min(nn,size(R_tilde,1));
                    assert(size(R_tilde,2) == nn);
                    kalman.steps{i}.R_tilde = R_tilde(1:ll,1:nn);
                    Q_T = Q';

                    o_c = [kalman.steps{i}.o ; kalman.steps{i_p_1}.c];
                    kalman.steps{i}.X = Q_T(1:ll,:) * [zeros(mm,size(D_i_p_1,2)) ; D_i_p_1];
                    kalman.steps{i}.o = Q_T(1:ll,:) * o_c;

                    kalman.steps{i_p_1}.D_tilde = Q_T(ll+1:end,:) * [zeros(mm,size(D_i_p_1,2)) ; D_i_p_1];
                    kalman.steps{i_p_1}.c = Q_T(ll+1:end,:) * o_c;
               end
            end

            for j=0:2:(l - 1)
                i = indices(j + 1);

                if (j == 0) % First index
                    kalman.steps{i}.R = kalman.steps{i}.R_tilde;
                    continue
                end
                R_tilde = kalman.steps{i}.R_tilde;
                B_i = kalman.steps{i}.B;
                D_i = kalman.steps{i}.D;

                nn = size(D_i,2);
                ll = size(D_i,1);

                [Q,R] = qr([D_i ; R_tilde]);
                assert( size(R,2) == nn );
                kalman.steps{i}.R = R(1:nn,1:nn);
                Q_T = Q';

                c_o = [kalman.steps{i}.c ; kalman.steps{i}.o];
                kalman.steps{i}.B_tilde = Q_T(1:nn,:)     * [B_i ; zeros(size(R_tilde,1),size(B_i,2))];
                kalman.steps{i}.c       = Q_T(1:nn,:)     * c_o;

                kalman.steps{i}.Z       = Q_T(nn+1:end,:) * [B_i ; zeros(size(R_tilde,1),size(B_i,2))];
                kalman.steps{i}.o       = Q_T(nn+1:end,:) * c_o;

                if (j + 1 ~= l) % Not last index
                    X = kalman.steps{i}.X;
                    kalman.steps{i}.Y       = Q_T(1:nn,:)     * [zeros(ll,size(X,2)) ; X];
                    kalman.steps{i}.X_tilde = Q_T(nn+1:end,:) * [zeros(ll,size(X,2)) ; X];
                end
            end

            for j=0:2:(l - 2)
                i = indices(j + 1);
                i_p_1 = indices(j + 2);

                D_tilde = kalman.steps{i_p_1}.D_tilde;
                C_i_p_1 = kalman.steps{i_p_1}.C;

                nn = size(D_tilde,2);

                [Q,R] = qr([D_tilde ; C_i_p_1]);
                assert(size(R,2)==nn);
                ll = nn;
                if size(R,1)<nn
                    ll = size(R,1); % this can happen in the first few steps if we do not have enough observations
                end
                kalman.steps{i_p_1}.C_tilde = R(1:ll,1:nn);
                Q_T = Q';

                %[i i_p_1 length(kalman.steps{i_p_1}.c) length(kalman.steps{i_p_1}.o)]
                %nn
                %ll
                %Q
                %R
                %Q_T
                c_o = [kalman.steps{i_p_1}.c ; kalman.steps{i_p_1}.o];
                %c_o
                kalman.steps{i_p_1}.c = Q_T(1:ll,:)     * c_o;
                kalman.steps{i_p_1}.o = Q_T(ll+1:end,:) * c_o;
                %if i_p_1==20; disp([111 length(kalman.steps{i_p_1}.o) length(kalman.steps{i_p_1}.c)]); end
                %kalman.steps{i_p_1}.c = Q_T(1:size(D_tilde,1),:) * c_o;
                %kalman.steps{i_p_1}.o = Q_T(size(D_tilde,1)+1:end,:) * c_o;

            end

            % fix the last index
            l_m_1 = indices(l - 1);
            last = indices(l);

            if (mod(l,2) == 1)
                C_tilde = kalman.steps{l_m_1}.C_tilde;
                o_i = kalman.steps{l_m_1}.c;

                Z = kalman.steps{last}.Z;
                c_i = kalman.steps{last}.o;


                [Q,R] = qr([C_tilde ; Z]);
                kalman.steps{l_m_1}.C_tilde = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                m = max(size(R,2),size(C_tilde,1)); % was just size(C_tilde,1), but this failed in test=[3 3 1 3]
                kalman.steps{l_m_1}.c = Q_T(1:m,:) * [o_i ; c_i];

                kalman.steps{last}.o = Q_T(m+1:end,:) * [o_i ; c_i];

                % if l_m_1==20 || last==20; disp([222 length(kalman.steps{last}.o) length(kalman.steps{l_m_1}.c)]); end
               % else
                %     C_tilde = kalman.steps{l_m_1}.C_tilde;
                %     o_i = kalman.steps{l_m_1}.c;
                %
                %     B_tilde = kalman.steps{last}.B_tilde;
                %     c_i = kalman.steps{last}.c;
                %
                %
                %     [Q,R] = qr([C_tilde ; B_tilde]);
                %     kalman.steps{l_m_1}.C_tilde = R(1:size(R,2),1:size(R,2));
                %     Q_T = Q';
                %
                %     kalman.steps{l_m_1}.c = Q_T(1:size(C_tilde,1),:) * [o_i ; c_i];
                %
                %     kalman.steps{last}.o = Q_T(size(C_tilde,1)+1:end,:) * [o_i ; c_i];
            end

            % Create new Kalman for the recursion
            % recursive_kalman = ParallelUltimateKalman();
            for j=0:2:(l - 2)
                i = indices(j + 1);
                i_p_1 = indices(j + 2);

                C_tilde = kalman.steps{i_p_1}.C_tilde;
                o = kalman.steps{i_p_1}.c;

                % D_tilde = kalman.steps{i + 1}.D_tilde;
                % C_tilde = kalman.steps{i + 1}.C_tilde;

                if (j ~= 0)
                    %disp('j != zero')
                    Z = kalman.steps{i}.Z;
                    X_tilde = kalman.steps{i}.X_tilde;
                    c = kalman.steps{i}.o;

                    % recursive_kalman.evolve(size(Z,2),X_tilde,-Z,c_i,CovarianceMatrix(ones(1,size(Z,1)),'w'))
                    kalman.steps{i_p_1}.B = Z;
                    kalman.steps{i_p_1}.D = X_tilde;

                    kalman.steps{i_p_1}.c = c;
                end

                % recursive_kalman.observe(C_tilde,o_i,CovarianceMatrix(ones(1,size(C_tilde,1)),'w'))
                kalman.steps{i_p_1}.C = C_tilde;
                kalman.steps{i_p_1}.o = o;
            end


            % The Recursion
            kalman.parallel_smooth(indices(2:2:l));

            % Now complete the recursive triangular solve & SelInv

            for j = 0:2:(l-1)
                i = indices(j + 1);

                R = kalman.steps{i}.R;

                if (j == 0)
                    i_p_1 = indices(j + 2);
                    x_i_p_1 = kalman.steps{i_p_1}.estimatedState;

                    X = kalman.steps{i}.X;
                    o_i = kalman.steps{i}.o;
                    kalman.steps{i}.estimatedState = R \ (o_i - X * x_i_p_1);

                    %fprintf('X0   set covariance for step %d; use %d\n',i,i_p_1);
                    S = eye(kalman.steps{i}.dimension) + X * kalman.steps{i_p_1}.S * X';
                    S = R \ ( S / (R') );
                    kalman.steps{i}.S = S;
                    kalman.steps{i}.S_off = - R \ ( X * kalman.steps{i_p_1}.S );
                    kalman.steps{i}.i_off = i_p_1;
                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix(S,'C');
                elseif (j ~= l - 1)
                    i_m_1 = indices(j);
                    i_p_1 = indices(j + 2);

                    x_i_p_1 = kalman.steps{i_p_1}.estimatedState;
                    x_i_m_1 = kalman.steps{i_m_1}.estimatedState;

                    B_tilde = kalman.steps{i}.B_tilde;
                    Y = kalman.steps{i}.Y;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - B_tilde * x_i_m_1 - Y * x_i_p_1);

                    %fprintf('X    set covariance for step %d; use %d, %d\n',i,i_p_1,i_m_1);
                    off_use = NaN;
                    Smm = kalman.steps{i_m_1}.S;
                    Spp = kalman.steps{i_p_1}.S;
                    if isfield(kalman.steps{i_p_1},'i_off') && kalman.steps{i_p_1}.i_off == i_m_1
                        off_use = i_p_1;
                        Smp = (kalman.steps{i_p_1}.S_off)';
                    end
                    if isfield(kalman.steps{i_p_1},'i1_off') && kalman.steps{i_p_1}.i1_off == i_m_1
                        off_use = -i_p_1; 
                        Smp = (kalman.steps{i_p_1}.S1_off)';
                    end
                    if isfield(kalman.steps{i_m_1},'i_off') && kalman.steps{i_m_1}.i_off == i_p_1 
                        off_use = i_m_1; 
                        Smp = (kalman.steps{i_m_1}.S_off);
                    end
                    if isfield(kalman.steps{i_m_1},'i1_off') && kalman.steps{i_m_1}.i1_off == i_p_1 
                        off_use = -i_m_1; 
                        Smp = (kalman.steps{i_m_1}.S1_off);
                    end
                    %fprintf('X    set covariance for step %d; use off %d\n',i,off_use);
                    %fprintf('B Y set covariance for step %d: %d %d (%d) \n',i, i_m_1, i_p_1,off_use);
                    S = eye(kalman.steps{i}.dimension) + (B_tilde * Smm * (B_tilde')) + (Y * Spp * (Y') + B_tilde * Smp * (Y') + Y * Smp' * (B_tilde'));
                    S = R \ ( S / (R') );

                    kalman.steps{i}.S = S;
                    kalman.steps{i}.i_off  = i_m_1;
                    kalman.steps{i}.S_off  = -R \ (B_tilde * Smm + Y * Smp' );
                    kalman.steps{i}.i1_off = i_p_1;
                    kalman.steps{i}.S1_off = -R \ (Y * Spp + B_tilde * Smp );

                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix(S,'C');
                else
                    i_m_1 = indices(j);

                    x_i_m_1 = kalman.steps{i_m_1}.estimatedState;

                    B_tilde = kalman.steps{i}.B_tilde;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - B_tilde * x_i_m_1);

                    %fprintf('Xl   set covariance for step %d; use %d\n',i,i_m_1);
                    S = eye(kalman.steps{i}.dimension) + B_tilde * kalman.steps{i_m_1}.S * B_tilde' ;
                    S = R \ ( S / (R') );
                    kalman.steps{i}.S = S;
                    kalman.steps{i}.S_off = - R \ B_tilde * kalman.steps{i_m_1}.S;
                    kalman.steps{i}.i_off = i_m_1;
                    kalman.steps{i}.estimatedCovariance = CovarianceMatrix(S,'C');
                end
            end
        end
    end

    methods (Access = public)
        function kalman = KalmanOddevenSmoother(options)
            if nargin<1
                options = struct();
            end
            kalman = kalman@KalmanExplicitRepresentation(options);
        end

        function smooth(kalman)
            for i = 1:length(kalman.steps)
                %i
                %kalman.steps{i}
                if kalman.steps{i}.step > 0
                    K_i = kalman.steps{i}.K;
                    kalman.steps{i}.B = - K_i.weigh(kalman.steps{i}.F);
                    kalman.steps{i}.c =   K_i.weigh(kalman.steps{i}.c);
                    kalman.steps{i}.D =   K_i.weigh(kalman.steps{i}.H);
                end
                if isfield(kalman.steps{i},'G')
                    C_i = kalman.steps{i}.C;
                    kalman.steps{i}.C = C_i.weigh(kalman.steps{i}.G);
                    kalman.steps{i}.o = C_i.weigh(kalman.steps{i}.o);
                else
                    kalman.steps{i}.C = [];
                    kalman.steps{i}.o = [];                    
                end

                n_i = kalman.steps{i}.dimension;
                kalman.steps{i}.estimatedCovariance = CovarianceMatrix( NaN * zeros(n_i,n_i), 'C' );
            end

            kalman.parallel_smooth();
        end


    end % methods
end