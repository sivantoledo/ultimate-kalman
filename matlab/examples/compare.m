function [fse,fce,sse,sce] = compare(kalman1, kalman2, seed, n, m, count, filter, smooth)
% COMPARE a function to compare the outputs of two Kalman filters/smoothers
%         mainly to test correctness of new implementations and variants.
%
%    COMPARE(kalman1, kalman2, seed, n, m, count, filter, smooth)
%      kalman1:    a just-constructed Kalman filter/smoother
%      kalman1:    a just-constructed Kalman filter/smoother
%      seed:       random-number generator seed
%      n:          dimension of state vectors
%      m:          dimension of observation vectors
%      count:      how many time steps to simulate and filter
%      filter:     should we compare filtered estimates and covariances
%      smooth:     should we compare smoothed estimates and covariances
%
% copyright 2024 Sivan Toledo

rng(1)

if ~isscalar(seed)
    test = seed;
    cov1type = n;
    cov2type = m;

    fse = NaN;
    fce = NaN;
    sse = NaN;
    sce = NaN;


    for i=1:size(test,1)
        n = test(i,1);
        l = test(i,2);
        m = test(i,3);
        f = test(i,4);
        r = test(i,5);

        if i==1
            kalman1.evolve(n);
            kalman2.evolve(n);
        else
            H = eye(l,n);
            [F,~] = qr( randn(max(l,n_i_m_1),max(l,n_i_m_1)) );
            F = F(1:l,1:n_i_m_1);
            c = randn(l,1);
            [K1,K2] = cov(l,cov1type,cov2type);
            %kalman1.evolve(n,H,F,c,CovarianceMatrix(Ke, 'C'))
            kalman1.evolve(n,H,F,c,K1)
            %kalman1.evolve(n,H,F,c,CovarianceMatrix(eye(l), 'C'))
            %kalman2.evolve(n,H,F,c,CovarianceMatrix(Ke, 'C'))
            kalman2.evolve(n,H,F,c,K2)
            %kalman2.evolve(n,H,F,c,CovarianceMatrix(eye(l), 'C'))
        end

        if m>0
            G = randn(m,n);
            o = randn(m,1);
            [C1,C2] = cov(m,cov1type,cov2type);
            %%Cf = randn(m,m);
            %%Ce = Cf' * Cf;
            %%L=chol(Ce,'lower');
            %%W=inv(L);
            %kalman1.observe(G,o,CovarianceMatrix(Ce, 'C'));
            kalman1.observe(G,o,C1);
            %kalman1.observe(G,o,CovarianceMatrix(eye(m), 'C'));
            %kalman2.observe(G,o,CovarianceMatrix(Ce, 'C'));
            kalman2.observe(G,o,C2);
            %kalman2.observe(G,o,CovarianceMatrix(eye(m), 'C'));
        else
            kalman1.observe([],[],CovarianceMatrix([],cov1type));
            kalman2.observe([],[],CovarianceMatrix([],cov1type));
        end

        n_i_m_1 = n;

        [u1,Z1] = kalman1.estimate();
        [u2,Z2] = kalman2.estimate();

        Z1 = Z1.explicit();
        Z2 = Z2.explicit();
        sz1 = sum(sum(Z1));
        sz2 = sum(sum(Z2));

        filteredStateErr(i) = max(abs((u1-u2)./u1));
        %filteredCovErr(1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
        filteredCovNorm       = max(max(abs(Z1)));
        filteredCovEltErr(i)   = max(max(abs((Z1-Z2)./Z1)));
        filteredCovNrmErr(i)   = max(max(abs((Z1-Z2)/filteredCovNorm)));

        if  isnan(sum(u1)) &&  isnan(sum(u2)); filteredStateErr(i)=   0; end
        if  isnan(sum(u1)) && ~isnan(sum(u2)); filteredStateErr(i)= inf; end
        if ~isnan(sum(u1)) &&  isnan(sum(u2)); filteredStateErr(i)= inf; end
        if  isnan(sz1) &&  isnan(sz2); filterCovNrmErr(i)=   0; end
        if  isnan(sz1) && ~isnan(sz2); filteredCovNrmErr(i)= inf; end
        if ~isnan(sz1) &&  isnan(sz2); filteredCovNrmErr(i)= inf; end

        %if filteredStateErr(i) > 1e-8
        %    i-1
        %    [u1 u2]
        %    error('large error');
        %end
    end

    %filteredStateErr
    fse = max(filteredStateErr);
    fce = max(filteredCovNrmErr);
    fprintf('max filtered relative error %.2e (elementwise) max filtered covariance error %.2e (normwise)\n', fse, fce);

    kalman1.smooth();
    kalman2.smooth();


    for i=1:size(test,1)

        [u1,Z1] = kalman1.estimate(i-1);
        [u2,Z2] = kalman2.estimate(i-1);

        Z1 = Z1.explicit();
        Z2 = Z2.explicit();
        sz1 = sum(sum(Z1));
        sz2 = sum(sum(Z2));

        %i
        %[Z1 Z2 Z1-Z2]

        smoothedStateErr(i) = max(abs((u1-u2)./u1));
        %filteredCovErr(1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
        smoothedCovNorm       = max(max(abs(Z1)));
        smoothedCovEltErr(i)   = max(max(abs((Z1-Z2)./Z1)));
        smoothedCovNrmErr(i)   = max(max(abs((Z1-Z2)/filteredCovNorm)));

        if  isnan(sum(u1)) &&  isnan(sum(u2)); smoothedStateErr(i)=   0; end
        if  isnan(sum(u1)) && ~isnan(sum(u2)); smoothedStateErr(i)= inf; end
        if ~isnan(sum(u1)) &&  isnan(sum(u2)); smoothedStateErr(i)= inf; end
        if  isnan(sz1) &&  isnan(sz2); smoothedCovNrmErr(i)=   0; end
        if  isnan(sz1) && ~isnan(sz2); smoothedCovNrmErr(i)= inf; end
        if ~isnan(sz1) &&  isnan(sz2); smoothedCovNrmErr(i)= inf; end
    end

    %smoothedCovNrmErr

    sse = max(smoothedStateErr);
    sce = max(smoothedCovNrmErr);
    fprintf('max smoothed relative error %.2e (elementwise) max smoothed covariance error %.2e (normwise)\n', sse, sce);

    return
end

rng(seed);

H = eye(n);
[F,~] = qr( randn(n,n), 0 );
%if m >= n
%    [G,~] = qr( randn(m,n), 0 );
%else
%    [G,~] = qr( randn(m,n) )
%end
G = randn(m,n);
c = zeros(n,1);
o = randn(m,1);

Kf = randn(n,n);
Ke = Kf' * Kf;
%K   = CovarianceMatrix(eye(n), 'C');
K   = CovarianceMatrix(Ke, 'C');
Cf = randn(m,m);
Ce = Cf' * Cf;
%C   = CovarianceMatrix(eye(m), 'C');
C   = CovarianceMatrix(Ce, 'C');

% cause lazy factorization
K.weigh(ones(n,1));
C.weigh(ones(m,1));

kalman1.evolve(n);
kalman2.evolve(n);

kalman1.observe(G,o,C);
kalman2.observe(G,o,C);

if filter
    [u1,Z1] = kalman1.estimate();
    [u2,Z2] = kalman2.estimate();
    %u1
    %u2
    %Z1e=Z1.explicit()
    %Ze2=Z2.explicit()

    filteredStateErr(1) = max(abs((u1-u2)./u1));
    %filteredCovErr(1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
    filteredCovNorm       = max(max(abs(Z1.explicit())));
    filteredCovEltErr(1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
    filteredCovNrmErr(1)   = max(max(abs((Z1.explicit()-Z2.explicit())/filteredCovNorm)));
end

for i=1:count

    if true
        [F,~] = qr( randn(n,n), 0 );
        G = randn(m,n);
        c = zeros(n,1);
        o = randn(m,1);
    end

    kalman1.evolve(n,H,F,c,K);
    kalman2.evolve(n,H,F,c,K);

    o = randn(m,1);
    kalman1.observe(G,o,C);
    kalman2.observe(G,o,C);

    if filter
        [u1,Z1] = kalman1.estimate();
        [u2,Z2] = kalman2.estimate();

        filteredStateErr(i+1) = max(abs((u1-u2)./u1));
        filteredCovNorm       = max(max(abs(Z1.explicit())));
        %filteredCovErr(i+1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
        filteredCovNorm       = max(max(abs(Z1.explicit())));
        filteredCovEltErr(i+1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
        filteredCovNrmErr(i+1)   = max(max(abs((Z1.explicit()-Z2.explicit())/filteredCovNorm)));
        %if filteredCovEltErr(i+1) > 1e-2
        %    fprintf("large error in covariance matrices of filtered estimates for step %d (%.2e)\n",i,filteredCovEltErr(i+1));
        %    [Z1.explicit() Z2.explicit()]
        %end
    end
end

if smooth
    kalman1.smooth();
    kalman2.smooth();

    for i=0:count
        [u1,Z1] = kalman1.estimate(i);
        [u2,Z2] = kalman2.estimate(i);

        smoothedStateErr(i+1) = max(abs((u1-u2)./u1));
        smoothedCovNorm       = max(max(abs(Z1.explicit())));
        smoothedCovEltErr(i+1)   = max(max(abs((Z1.explicit()-Z2.explicit())./Z1.explicit())));
        smoothedCovNrmErr(i+1)   = max(max(abs((Z1.explicit()-Z2.explicit())/smoothedCovNorm)));

        %if smoothedCovEltErr(i+1) > 1e-2
        %    fprintf("large error in covariance matrices of smoothed estimates for step %d (%.2e)\n",i,smoothedCovEltErr(i+1));
        %    [Z1.explicit() Z2.explicit()]
        %end
        %[Z1.explicit() Z2.explicit()]
        %if smoothedCovNrmErr(i+1) > 1e-12
        %    fprintf('err in %d (%.2e)\n',i+1,smoothedCovNrmErr(i+1))
        %z1e=Z1.explicit();
        %z2e=Z2.explicit();
        %[z1e z2e z1e-z2e]
        %end
    end
end

close all;
figure

if filter
    subplot(2,2,1);
    plot(filteredStateErr);
    title('Fitered State Elementwise Errors');

    subplot(2,2,2);
    semilogy(filteredCovEltErr);
    hold on;
    semilogy(filteredCovNrmErr);
    hold off
    legend('Elementwise','Normwise')
    title('Filtered Covariance Errors');

    fprintf('max filtered relative error %.2e (elementwise) max filtered covariance error %.2e (elementwise) %.2e (normwise)\n' ...
        ,max(filteredStateErr), max(filteredCovEltErr), max(filteredCovNrmErr));
end

if smooth
    subplot(2,2,3);
    plot(smoothedStateErr);
    title('Smoothed State Elementwise Errors');

    subplot(2,2,4);
    semilogy(smoothedCovEltErr);
    hold on;
    semilogy(smoothedCovNrmErr);
    hold off
    legend('Elementwise','Normwise')
    title('Smoothed Covariance Errors');

    fprintf('max smoothed relative error %.2e (elementwise) max filtered covariance error %.2e (elementwise) %.2e (normwise)\n' ...
        ,max(smoothedStateErr), max(smoothedCovEltErr), max(smoothedCovNrmErr));
end

    function [cov1,cov2] = cov(d,cov1type,cov2type)
        if cov1type=='w' || cov2type=='w'
            covdiag = true;
        else
            covdiag = false;
        end

        clear cov1 cov2;
        if covdiag
            w   = randn(d,1);
            cov = diag(w.^-2);
            W   = diag(w);
            if cov1type=='W'
                cov1=CovarianceMatrix(W, 'W');
            end
            if cov1type=='w'
                cov1=CovarianceMatrix(w, 'w');
            end
            if cov1type=='C'
                cov1=CovarianceMatrix(cov, 'C');
            end
            if cov2type=='W'
                cov2=CovarianceMatrix(W, 'W');
            end
            if cov2type=='w'
                cov2=CovarianceMatrix(w, 'w');
            end
            if cov2type=='C'
                cov2=CovarianceMatrix(cov, 'C');
            end
        else
            fact = randn(d,d);
            cov  = fact' * fact;
            L=chol(cov,'lower');
            W=inv(L);
            if cov1type=='W'
                cov1=CovarianceMatrix(W, 'W');
            end
            if cov1type=='C'
                cov1=CovarianceMatrix(cov, 'C');
            end
            if cov2type=='W'
                cov2=CovarianceMatrix(W, 'W');
            end
            if cov2type=='C'
                cov2=CovarianceMatrix(cov, 'C');
            end
        end
    end
end