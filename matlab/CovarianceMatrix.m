classdef CovarianceMatrix < handle
    % CovarianceMatrix   Representation of a covariance matric C for use in 
    %                    generalized least squares problems.
    %                    We denote W'*W = inv(C)
    %
    % CovarianceMatrix methods:
    %    weigh    - Multiplies a given matrix or column vector by W
    %    explicit - Returns an explicit dense covariance matrix C 
    %
    % Constructors:
    %    CovarianceMatrix(Z,type) where Z is a representation of C and type
    %                             defines how Z is related to C.
    %
    % type=='C' (covariance):           Z=C,         L=chol(C,'lower'), W=inv(L)
    % type=='I' (inverse covariance):   Z=inv(C),    L=chol(C,'lower')  W=L'
    % type=='M' (sqrt of covariance):   Z*Z'=C, e.g. Z=chol(C,'lower')  W=inv(Z)
    % type=='W' (inverse factor):       Z'*Z=inv(C),                    W=Z
    % type=='w' (diagonal inv. factor): diag(Z)*diag(Z)=inv(C),                    W=Z
    %
    % Copyright 2020-2024 Sivan Toledo, Tel Aviv University.

 
    properties
        type; 
        invW;
        W;
        w;
        C; % for lazy factorization
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
                    % Sivan July 2024, lazy factorization for conventional
                    % filters.
                    %cov.type = 'U';
                    %cov.invW = chol(Z);
                    cov.type = 'C';
                    cov.C    = Z;
                case 'I'       % the inverse of the covariance matrix 
                    cov.type = 'W';
                    cov.W = (chol(Z,'lower'))';
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
                case {'U', 'M', 'F'}       
                    WA = (cov.invW)\A;
                case 'W'
                    WA = (cov.W)*A;
                case 'w'
                    WA = diag(cov.w)*A;
                case 'C'
                    cov.type = 'F'; % lowe-Cholesky factored
                    L = chol(cov.C,'lower'); 
                    cov.invW = L; 
                    WA = (cov.invW)\A;
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
                case 'C'
                    R = cov.C;
                case 'F'
                    % return just one, for the mex interface
                    %R = { cov.C, cov.invW };    
                    R = cov.invW;    
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
                case {'C', 'F'}
                    C = cov.C;
                otherwise
                    error('illegal covariance type');
            end
        end

    end
end
