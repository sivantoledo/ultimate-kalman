function t = compare(kalman1, kalman2, seed, n, m, count, filter, smooth)
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

rng(seed);

H = eye(n);
[F,~] = qr( randn(n,n), 0 );
[G,~] = qr( randn(m,n), 0 );
c = zeros(n,1);
o = randn(m,1);

K   = CovarianceMatrix(eye(n), 'C');
C   = CovarianceMatrix(eye(m), 'C');

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

        %if smoothedCovNrmErr(i+1) > 1e-12
        %    z1e=Z1.explicit()
        %    z2e=Z2.explicit()
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