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

            if (length(indices) == 1)
                singleStep = kalman.steps{indices(1)};
                singleStep.estimatedState = singleStep.C \ singleStep.o;

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

                    [Q,R_tilde] = qr([C_i ; B_i_p_1]);
                    kalman.steps{i}.R_tilde = R_tilde(1:size(R_tilde,2),1:size(R_tilde,2));
                    Q_T = Q';

                    o_c = [kalman.steps{i}.o ; kalman.steps{i_p_1}.c];
                    kalman.steps{i}.X = Q_T(1:size(C_i,1),:) * [zeros(size(C_i)) ; D_i_p_1];
                    kalman.steps{i}.o = Q_T(1:size(C_i,1),:) * o_c;

                    kalman.steps{i_p_1}.D_tilde = Q_T(size(C_i,1)+1:end,:) * [zeros(size(C_i)) ; D_i_p_1];
                    kalman.steps{i_p_1}.c = Q_T(size(C_i,1)+1:end,:) * o_c;
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

                [Q,R] = qr([D_i ; R_tilde]);
                kalman.steps{i}.R = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                c_o = [kalman.steps{i}.c ; kalman.steps{i}.o];
                kalman.steps{i}.B_tilde = Q_T(1:size(B_i,1),:) * [B_i ; zeros(size(R_tilde))];
                kalman.steps{i}.c = Q_T(1:size(B_i,1),:) * c_o;

                kalman.steps{i}.Z = Q_T(size(B_i,1)+1:end,:) * [B_i ; zeros(size(R_tilde))];
                kalman.steps{i}.o = Q_T(size(B_i,1)+1:end,:) * c_o;

                if (j + 1 ~= l) % Not last index
                    X = kalman.steps{i}.X;
                    kalman.steps{i}.Y = Q_T(1:size(B_i,1),:) * [zeros(size(B_i)) ; X];

                    kalman.steps{i}.X_tilde = Q_T(size(B_i,1)+1:end,:) * [zeros(size(B_i)) ; X];
                end
            end

            for j=0:2:(l - 2)
                i = indices(j + 1);
                i_p_1 = indices(j + 2);

                D_tilde = kalman.steps{i_p_1}.D_tilde;
                C_i_p_1 = kalman.steps{i_p_1}.C;

                [Q,R] = qr([D_tilde ; C_i_p_1]);
                kalman.steps{i_p_1}.C_tilde = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                c_o = [kalman.steps{i_p_1}.c ; kalman.steps{i_p_1}.o];
                kalman.steps{i_p_1}.c = Q_T(1:size(D_tilde,1),:) * c_o;

                kalman.steps{i_p_1}.o = Q_T(size(D_tilde,1)+1:end,:) * c_o;

            end


            %fix the last index
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

                kalman.steps{l_m_1}.c = Q_T(1:size(C_tilde,1),:) * [o_i ; c_i];

                kalman.steps{last}.o = Q_T(size(C_tilde,1)+1:end,:) * [o_i ; c_i];
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

            for j = 0:2:(l-1)
                i = indices(j + 1);


                R = kalman.steps{i}.R;

                if (j == 0)
                    i_p_1 = indices(j + 2);
                    x_i_p_1 = kalman.steps{i_p_1}.estimatedState;

                    X = kalman.steps{i}.X;
                    o_i = kalman.steps{i}.o;
                    kalman.steps{i}.estimatedState = R \ (o_i - X * x_i_p_1);
                elseif (j ~= l - 1)
                    i_m_1 = indices(j);
                    i_p_1 = indices(j + 2);

                    x_i_p_1 = kalman.steps{i_p_1}.estimatedState;
                    x_i_m_1 = kalman.steps{i_m_1}.estimatedState;

                    B_tilde = kalman.steps{i}.B_tilde;
                    Y = kalman.steps{i}.Y;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - B_tilde * x_i_m_1 - Y * x_i_p_1);
                else
                    i_m_1 = indices(j);

                    x_i_m_1 = kalman.steps{i_m_1}.estimatedState;

                    B_tilde = kalman.steps{i}.B_tilde;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - B_tilde * x_i_m_1);
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
                end

                n_i = kalman.steps{i}.dimension;
                kalman.steps{i}.estimatedCovariance = CovarianceMatrix( NaN * zeros(n_i,n_i), 'C' );
            end

            kalman.parallel_smooth();
        end


    end % methods
end