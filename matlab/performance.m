function t = performance(kalman_factory, seed, sizes, count, decimation)

if nargin < 1
    seed = 1;
end

rng(seed);

%decimation = min(1000,count);

%t = NaN*zeros(ceil(count/100),length(sizes));

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

    kalman = kalman_factory();

    if true
        t(:,s) = kalman.perftest([],F,c,K,G,o,C,count,decimation);
    else
        j = 1;

        tic;

        for i=1:count
            if rem(i,100) == 1
                j=j+1;
                t(j,s) = toc;
            end

            kalman.evolve(n,[],F,c,K);
            kalman.observe(G,randn(m,1),C);
            filtered = kalman.estimate();
            if filtered(1)==7 % just to make sure we use the result, to avoid deleting it as dead code
                disp('dummy');
            end
            kalman.forget();
        end
    end
end

%t = t(10:end,:);

%r = (t(2:end,:) - t(1:end-1,:))/100;
r = 1000*t; % r*1000; % milliseconds

%r

medians = median(r)
size(r)
sizes

close all
figure
axis square
set(gca,'Box','on');;
hold on;
%plot(r,'k-','LineWidth',1);
plot(decimation*(0:size(r,1)-1),r);
xlim([0 decimation*(size(r,1)-1)]);
ylim([0 1.33*max(medians)])
legend(num2str(sizes'),'Location','SouthEast');
xlabel('step');
ylabel('time (ms)')
hold off;
%exportgraphics(gca,'../outputs/stress_6_96.pdf');



