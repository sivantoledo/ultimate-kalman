classdef CovarianceMatrix < handle

    % covType=='C' (covariance):         L*L'=Z, W=inv(L)
    % covType=='I' (inverse covariance): L*L'=Z, W=L'
    % covType=='M' (sqrt of covariance):         W=inv(Z)
    % covType=='W' (inverse sqrt):               W=Z
 
    properties
        type; 
        invW;
        W;
        w;
    end

    methods 
        function cov = CovarianceMatrix(Z,type)
            %CovarianceMatrix constructor
            %
            % cov = CovarianceMatrix(Z,type)
            %    Z is a real matrix or vector
            %    type describes how Z represents a covariance matrix C
            %    'C': C = Z
            %    'I': C = inv(Z)
            %    'W': C = inv( Z' * Z )
            %    'w': C = diag( Z.^-2 ) = inv( diag(Z)' * diag(Z) )
            %    'M': C = Z' * Z
            %    'U': C = Z * Z' and Z is upper triangular
            %
            cov = cov@handle();
            switch type
                case 'C'       % an explicit covariance matrix
                    cov.type = 'U';
                    cov.invW = chol(Z);
                case 'I'       % the inverse of the covariance matrix 
                    cov.type = 'W';
                    cov.W = (chol(Z))';
                case 'M'       % a factor (perhaps Cholesky factor) of the covariance matrix 
                    cov.type = 'M'; % hopefully never happens...
                    if (istriu(Z)); cov.type = 'U'; end;
                    cov.invW = Z;
                case 'W'       % an inverse factor W such that W'*W=C
                    cov.type = 'W';
                    cov.W = Z;
                case 'w'       % an inverse factor W such that W'*W=C
                    cov.type = 'w';
                    cov.w = Z;
                otherwise
                    type
                    error('illegal covariance-matrix representation');
            end
        end

        function WA = weigh(cov,A)
            switch cov.type
                case {'U', 'M'}       
                   WA = (cov.invW)\A;
                case 'W'
                    WA = (cov.W)*A;
                case 'w'
                    WA = diag(cov.w)*A;
                otherwise
                    error('illegal covariance type');
            end
        end

        function [R,type] = rep(cov)
            type = cov.type;
            switch cov.type
                case {'U', 'M'}
                    R = cov.invW;
                case 'W'
                    R = cov.W;
                case 'w'
                    R = cov.w;
                otherwise
                    error('illegal covariance type');
            end
        end

        function C = explicit(cov)
            switch cov.type
                case {'U', 'M'}
                    C = (cov.invW)' * (cov.invW);
                case 'W'
                    C = inv( (cov.W)' * (cov.W) );
                case 'w'
                    C = diag(cov.w.^-2);
                otherwise
                    error('illegal covariance type');
            end
        end

    end
end
