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
                singleStep.estimatedState = singleStep.G \ singleStep.o;

                kalman.steps{indices(1)} = singleStep;
                return
            end

            l = length(indices);

            for j=0:2:(l - 1)
                i = indices(j + 1);

                if (j + 1 == l) % Last index
                    G_i = kalman.steps{i}.G;
                    [Q,R_tilde] = qr(G_i);
                    kalman.steps{i}.R_tilde = R_tilde;

                    kalman.steps{i}.o = Q' * kalman.steps{i}.o;
                else
                    i_p_1 = indices(j + 2);
                    G_i = kalman.steps{i}.G;
                    F_i_p_1 = kalman.steps{i_p_1}.F;
                    H_i_p_1 = kalman.steps{i_p_1}.H;

                    [Q,R_tilde] = qr([G_i ; F_i_p_1]);
                    kalman.steps{i}.R_tilde = R_tilde(1:size(R_tilde,2),1:size(R_tilde,2));
                    Q_T = Q';

                    o_c = [kalman.steps{i}.o ; kalman.steps{i_p_1}.c];
                    kalman.steps{i}.X = Q_T(1:size(G_i,1),:) * [zeros(size(G_i)) ; H_i_p_1];
                    kalman.steps{i}.o = Q_T(1:size(G_i,1),:) * o_c;

                    kalman.steps{i_p_1}.H_tilde = Q_T(size(G_i,1)+1:end,:) * [zeros(size(G_i)) ; H_i_p_1];
                    kalman.steps{i_p_1}.c = Q_T(size(G_i,1)+1:end,:) * o_c;
                end
            end

            for j=0:2:(l - 1)
                i = indices(j + 1);

                if (j == 0) % First index
                    kalman.steps{i}.R = kalman.steps{i}.R_tilde;
                    continue
                end
                R_tilde = kalman.steps{i}.R_tilde;
                F_i = kalman.steps{i}.F;
                H_i = kalman.steps{i}.H;

                [Q,R] = qr([H_i ; R_tilde]);
                kalman.steps{i}.R = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                c_o = [kalman.steps{i}.c ; kalman.steps{i}.o];
                kalman.steps{i}.F_tilde = Q_T(1:size(F_i,1),:) * [F_i ; zeros(size(R_tilde))];
                kalman.steps{i}.c = Q_T(1:size(F_i,1),:) * c_o;

                kalman.steps{i}.Z = Q_T(size(F_i,1)+1:end,:) * [F_i ; zeros(size(R_tilde))];
                kalman.steps{i}.o = Q_T(size(F_i,1)+1:end,:) * c_o;

                if (j + 1 ~= l) % Not last index
                    X = kalman.steps{i}.X;
                    kalman.steps{i}.Y = Q_T(1:size(F_i,1),:) * [zeros(size(F_i)) ; X];

                    kalman.steps{i}.X_tilde = Q_T(size(F_i,1)+1:end,:) * [zeros(size(F_i)) ; X];
                end
            end

            for j=0:2:(l - 2)
                i = indices(j + 1);
                i_p_1 = indices(j + 2);

                H_tilde = kalman.steps{i_p_1}.H_tilde;
                G_i_p_1 = kalman.steps{i_p_1}.G;

                [Q,R] = qr([H_tilde ; G_i_p_1]);
                kalman.steps{i_p_1}.G_tilde = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                c_o = [kalman.steps{i_p_1}.c ; kalman.steps{i_p_1}.o];
                kalman.steps{i_p_1}.c = Q_T(1:size(H_tilde,1),:) * c_o;

                kalman.steps{i_p_1}.o = Q_T(size(H_tilde,1)+1:end,:) * c_o;

            end


            %fix the last index
            l_m_1 = indices(l - 1);
            last = indices(l);

            if (mod(l,2) == 1)
                G_tilde = kalman.steps{l_m_1}.G_tilde;
                o_i = kalman.steps{l_m_1}.c;

                Z = kalman.steps{last}.Z;
                c_i = kalman.steps{last}.o;


                [Q,R] = qr([G_tilde ; Z]);
                kalman.steps{l_m_1}.G_tilde = R(1:size(R,2),1:size(R,2));
                Q_T = Q';

                kalman.steps{l_m_1}.c = Q_T(1:size(G_tilde,1),:) * [o_i ; c_i];

                kalman.steps{last}.o = Q_T(size(G_tilde,1)+1:end,:) * [o_i ; c_i];
                % else
                %     G_tilde = kalman.steps{l_m_1}.G_tilde;
                %     o_i = kalman.steps{l_m_1}.c;
                %
                %     F_tilde = kalman.steps{last}.F_tilde;
                %     c_i = kalman.steps{last}.c;
                %
                %
                %     [Q,R] = qr([G_tilde ; F_tilde]);
                %     kalman.steps{l_m_1}.G_tilde = R(1:size(R,2),1:size(R,2));
                %     Q_T = Q';
                %
                %     kalman.steps{l_m_1}.c = Q_T(1:size(G_tilde,1),:) * [o_i ; c_i];
                %
                %     kalman.steps{last}.o = Q_T(size(G_tilde,1)+1:end,:) * [o_i ; c_i];
            end



            % Create new Kalman for the recursion
            % recursive_kalman = ParallelUltimateKalman();
            for j=0:2:(l - 2)
                i = indices(j + 1);
                i_p_1 = indices(j + 2);

                G_tilde = kalman.steps{i_p_1}.G_tilde;
                o = kalman.steps{i_p_1}.c;

                % H_tilde = kalman.steps{i + 1}.H_tilde;
                % G_tilde = kalman.steps{i + 1}.G_tilde;

                if (j ~= 0)
                    Z = kalman.steps{i}.Z;
                    X_tilde = kalman.steps{i}.X_tilde;
                    c = kalman.steps{i}.o;

                    % recursive_kalman.evolve(size(Z,2),X_tilde,-Z,c_i,CovarianceMatrix(ones(1,size(Z,1)),'w'))
                    kalman.steps{i_p_1}.F = Z;
                    kalman.steps{i_p_1}.H = X_tilde;

                    kalman.steps{i_p_1}.c = c;
                end

                % recursive_kalman.observe(G_tilde,o_i,CovarianceMatrix(ones(1,size(G_tilde,1)),'w'))
                kalman.steps{i_p_1}.G = G_tilde;
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

                    F_tilde = kalman.steps{i}.F_tilde;
                    Y = kalman.steps{i}.Y;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - F_tilde * x_i_m_1 - Y * x_i_p_1);
                else
                    i_m_1 = indices(j);

                    x_i_m_1 = kalman.steps{i_m_1}.estimatedState;

                    F_tilde = kalman.steps{i}.F_tilde;
                    c = kalman.steps{i}.c;

                    kalman.steps{i}.estimatedState = R \ (c - F_tilde * x_i_m_1);
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
                    kalman.steps{i}.F = - K_i.weigh(kalman.steps{i}.F);
                    kalman.steps{i}.c =   K_i.weigh(kalman.steps{i}.c);
                    kalman.steps{i}.H =   K_i.weigh(kalman.steps{i}.H);
                end
                if isfield(kalman.steps{i},'G')
                    C_i = kalman.steps{i}.C;
                    kalman.steps{i}.G = C_i.weigh(kalman.steps{i}.G);
                    kalman.steps{i}.o = C_i.weigh(kalman.steps{i}.o);
                end

                n_i = kalman.steps{i}.dimension;
                kalman.steps{i}.estimatedCovariance = CovarianceMatrix( NaN * zeros(n_i,n_i), 'C' );
            end

            kalman.parallel_smooth();
        end


    end % methods
end