classdef KalmanFilterSmoother < KalmanBase
    % KalmanFilterSmoother   A conventional Kalman filter and RTS (Rauch,
    %                        Tung, and Striebel) smoother. The
    %                        implementation works both as a linear
    %                        filter/smoother or an extended one, depending
    %                        on whether the evolution and observation
    %                        operators are matrices or functions.
    %
    % KalmanFilterSmoother inherits from KalmanBase.
    %
    % KalmanFilterSmoother fields (properties): none
    %  
    % KalmanFilterSmoother implements these abstract methods:
    %    evolve   - Evolve the state using a linear matrix equation 
    %    observe  - Provide observations of the current state 
    %    smooth   - Compute smooth estimates of all the stored states
    % 
    % Copyright 2024 Sivan Toledo, Tel Aviv University.

    properties (Access = protected)
    end

    methods (Access = protected)
        function rollbackStepWithIndex(kalman,i)
            kalman.current           = [];
            kalman.current.dimension = kalman.steps{ i }.dimension;
            kalman.current.step      = kalman.steps{ i }.step;
            kalman.current.predictedState      = kalman.steps{ i }.predictedState;
            kalman.current.predictedCovariance = kalman.steps{ i }.predictedCovariance;
            kalman.current.F                   = kalman.steps{ i }.F;
            kalman.current.H                   = kalman.steps{ i }.H;
        end
    end

    methods (Access = public)
        function kalman = KalmanFilterSmoother(options)
            if nargin==0
                options = struct();
            end
            kalman = kalman@KalmanBase(options);
        end

        function evolve(kalman,n_i,H_i,F_i,c_i,K_i)
            %EVOLVE   Evolve the state using the given linear recurrence.
            %   kalman.EVOLVE(H_i,F_i,c_i,K_i) evolves the state using the recurrence
            %                H_i * u_i = F_i * u_{i-1} + c_i + epsilon
            %   where C_i is the covariance matrix of the error term epsilon.
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

            %fprintf('evolve %d\n',length(kalman.steps))

            kalman.current = []; % clear
            kalman.current.dimension = n_i;
            if isempty(kalman.steps) % starting the first step
                kalman.current.step = 0;
                return
            end

                
            ptr_imo = length(kalman.steps);                                % pointer to row i-1 in the cell array
            kalman.current.step  = kalman.steps{ptr_imo}.step + 1;

            if isempty(H_i)
                H_i = eye(size(F_i,1));
            end
            %nargin
            %isempty(H_i)
            %norm(H_i-eye(size(H_i,1)>0))
            
            %if ~isempty(H_i) && (size(H_i,1) ~= size(H_i,2) || norm(H_i-eye(size(H_i,1)))>0)
            %   error('This implementation of UltimateKalman only accepts H_i==[] or H_i=eye');
            %end

            % keep F_i for smoothing
            kalman.current.F                   = F_i;
            kalman.current.H                   = H_i;

            % predictions
            if isa(F_i,'function_handle')
                [predictedState,jacobian] = F_i( kalman.steps{ptr_imo}.assimilatedState );
            else
                predictedState = F_i * kalman.steps{ptr_imo}.assimilatedState;
                jacobian       = F_i;
            end
            


            %kalman.current.predictedState      = predictedState + c_i;
            %kalman.current.predictedCovariance = jacobian * kalman.steps{ptr_imo}.assimilatedCovariance * jacobian' + K_i.explicit(); 
            Hinv = inv(H_i);
            kalman.current.predictedState      = Hinv * (predictedState + c_i);
            kalman.current.predictedCovariance = Hinv * jacobian * kalman.steps{ptr_imo}.assimilatedCovariance * jacobian' * Hinv' + Hinv * K_i.explicit() * Hinv'; 

            kalman.current.estimatedState      = kalman.current.predictedState;
            kalman.current.estimatedCovariance = CovarianceMatrix( kalman.current.predictedCovariance, 'C');
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

            if kalman.current.step == 0
                
                if isa(G_i,'function_handle')
                    error('in an extended Kalman filter, state 0 must be observed linearly')
                end

                W_i_G_i = C_i.weigh(G_i);
                W_i_o_i = C_i.weigh(o_i);
                [Q,R] = qr(W_i_G_i,0); 
                kalman.current.assimilatedState      = R \ ( Q' * W_i_o_i );
                kalman.current.assimilatedCovariance = inv(R'*R); 

                kalman.current.estimatedState      = kalman.current.assimilatedState;
                kalman.current.estimatedCovariance = CovarianceMatrix( kalman.current.assimilatedCovariance, 'C');
                
                kalman.steps{ length(kalman.steps)+1 } = kalman.current;
                return;
            end

            if nargin<4 || isempty(o_i) % no observations, pass []
                %fprintf('observe step %d, no observations\n',length(kalman.steps))
                kalman.current.assimilatedState      = kalman.current.predictedState;
                kalman.current.assimilatedCovariance = kalman.current.predictedCovariance; 
            else % there are observations
                if isa(G_i,'function_handle')
                    [predictedObservations,jacobian] = G_i(kalman.current.predictedState);
                else
                    predictedObservations = G_i * kalman.current.predictedState;
                    jacobian = G_i;
                end
                S          = jacobian * kalman.current.predictedCovariance * jacobian' + C_i.explicit;
                gain       = kalman.current.predictedCovariance * jacobian' / S;
                innovation = o_i - predictedObservations;

                % the formula in Sarkka's book for the assimilated
                % covariance is 
                % predictedCovariance - ...
                %     gain * ( G_i * predictedCovariance * G_i' + C_i ) * gain' 
                % which simplifies to the expression here, which is the one
                % shown in Strang & Borre
                kalman.current.assimilatedState      = kalman.current.predictedState + gain * innovation;
                kalman.current.assimilatedCovariance = kalman.current.predictedCovariance - gain * jacobian * kalman.current.predictedCovariance; 
                %kalman.current.assimilatedCovariance = kalman.current.predictedCovariance - gain * S * gain'; 
            end % of we have observations

            kalman.current.estimatedState      = kalman.current.assimilatedState;
            kalman.current.estimatedCovariance = CovarianceMatrix( kalman.current.assimilatedCovariance, 'C');

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
            
            % Warning('smoothed covariances seem wrong');

            l = length(kalman.steps);

            kalman.steps{l}.smoothedState      = kalman.steps{l}.assimilatedState;
            kalman.steps{l}.smoothedCovariance = kalman.steps{l}.assimilatedCovariance;

            for i=l-1:-1:1
                nextSmoothEstimate      = kalman.steps{i+1}.smoothedState;
                nextPredictedEstimate   = kalman.steps{i+1}.predictedState;
                nextPredictedCovariance = kalman.steps{i+1}.predictedCovariance;
                nextSmoothedCovariance  = kalman.steps{i+1}.smoothedCovariance;
                F_ipo = kalman.steps{i+1}.F;
                if isa(F_ipo,'function_handle')
                    % this exact jacobian was already computed in evolve,
                    % makes more sense to keep it
                    [~,nextEvolutionMatrix] = F_ipo( kalman.steps{i}.assimilatedState ) ;
                else
                    nextEvolutionMatrix     = F_ipo;
                end
                H_ipo = kalman.steps{i+1}.H;
                nextEvolutionMatrix = H_ipo \ nextEvolutionMatrix ;

                assimilatedCovariance = kalman.steps{i}.assimilatedCovariance;
                assimilatedEstimate   = kalman.steps{i}.assimilatedState;

                %backwardInnovation = assimilatedCovariance * nextEvolutionMatrix' * inv(nextPredictedCovariance);
                backwardInnovation = assimilatedCovariance * nextEvolutionMatrix' / nextPredictedCovariance;

                % formulas for the RTS smoother and its covariance matrices
                % from Sarkka's book.
                kalman.steps{i}.smoothedState = assimilatedEstimate + backwardInnovation * ( nextSmoothEstimate - nextPredictedEstimate );
                kalman.steps{i}.smoothedCovariance = assimilatedCovariance ... 
                    + backwardInnovation * ( ...
                                nextSmoothedCovariance - nextPredictedCovariance ...
                            ) * backwardInnovation';

                kalman.steps{i}.estimatedState      = kalman.steps{i}.smoothedState;
                kalman.steps{i}.estimatedCovariance = CovarianceMatrix( kalman.steps{i}.smoothedCovariance, 'C' );
            end

        end
    end % methods
end