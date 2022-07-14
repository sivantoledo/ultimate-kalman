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
        end

        function WA = weigh(cov,A)
            if cov.type=='C'       
                WA = (cov.invW)\A;
            end
            if cov.type=='I'
                WA = (cov.Wtrans)'*A;
            end
            if cov.type=='M'
                invW = Z;
                WA = (cov.invW)\A;
            end
            if cov.type=='W'
                WA = (cov.W)*A;
            end
        end
    end
end
