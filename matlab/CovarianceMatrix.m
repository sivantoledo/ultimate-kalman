classdef CovarianceMatrix < handle

    % covType=='C' (covariance):         L*L'=Z, W=inv(L)
    % covType=='I' (inverse covariance): L*L'=Z, W=L'
    % covType=='M' (sqrt of covariance):         W=inv(Z)
    % covType=='W' (inverse sqrt):               W=Z
 
    properties
        type; 
        invW;
        Wtrans;
        W;
        w;
    end

    methods 
        function cov = CovarianceMatrix(Z,type)
            cov = cov@handle();
            cov.type = type;
            if type=='C'       % an explicit covariance matrix
                cov.invW = chol(Z);
            end
            if type=='I'       % the inverse of the covariance matrix 
                cov.Wtrans = chol(Z);
            end
            if type=='M'       % a factor (perhaps Cholesky factor) of the covariance matrix 
                cov.invW = Z;
            end
            if type=='W'       % an inverse factor W such that W'*W=C
                cov.W = Z;
            end
            if type=='w'       % an inverse factor W such that W'*W=C
                cov.w = Z;
            end
        end

        function WA = weigh(cov,A)
            if cov.type=='C'       
                WA = (cov.invW)\A;
            end
            if cov.type=='I'
                WA = (cov.Wtrans)'*A;
            end
            if cov.type=='M'
                WA = (cov.invW)\A;
            end
            if cov.type=='W'
                WA = (cov.W)*A;
            end
            if cov.type=='w'
                WA = diag(cov.w)*A;
            end
        end

        function [R,type] = rep(cov)
            if cov.type=='C'       % an explicit covariance matrix
                type = 'U'; % an upper triangular matrix, WA = R \ A;
                R = cov.invW;
            end
            if cov.type=='I'       % the inverse of the covariance matrix
                type = 'W';
                R = (cov.Wtrans)';
            end
            if cov.type=='M'       % a factor (perhaps Cholesky factor) of the covariance matrix 
                R = cov.invW;
                type = '?';
                if (istriu(R)); type = 'U'; end;
                if (istriu(l)); type = 'L'; end;
            end
            if cov.type=='W'       % an inverse factor W such that W'*W=C
                type = 'W';
                R = cov.W;
            end
            if cov.type=='w'       % an inverse factor W such that W'*W=C
                type = 'w';
                R = cov.w;
            end
        end

        function C = explicit(cov)
            if cov.type=='C'       
                C = (cov.invW)' * (cov.invW);
            end
            if cov.type=='I'
                C = inv(cov.Wtrans * cov.Wtrans');
            end
            if cov.type=='M'
                C = (cov.invW)' * (cov.invW);
            end
            if cov.type=='W'
                C = inv(cov.W'  *cov.W);
            end
            if cov.type=='w'
                C = inv(diag(cov.w.^2));
            end
        end

    end
end
