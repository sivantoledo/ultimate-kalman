function t = performance_smooth(seed, sizes, count)

if nargin < 1
    seed = 1;
end

rng(seed);

q = 1;

for s = 1:length(sizes)

    n = sizes(s)

    l = n;
    m = n;

    [F,~] = qr( randn(n,n) );
    [G,~] = qr( randn(m,n) );
    c = zeros(n,1);
    o = randn(m,1);

    K   = CovarianceMatrix(eye(n), 'C');
    C   = CovarianceMatrix(eye(m), 'C');

    kalman = UltimateKalmanNative;

    tstart=tic;
    t = kalman.perftest_smooth([],F,c,K,G,o,C,count);
    tend=toc(tstart);

    [n count tend tend/count t]
end

