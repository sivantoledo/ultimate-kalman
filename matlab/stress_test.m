function t = stress_test(seed, sizes , k)

if nargin < 1
    seed = 1;
end

rng(seed);

t = zeros(ceil(k/100),length(sizes));

for s = 1:length(sizes)

    n = sizes(s)



    l = n;
    m = n;

    [F,~] = qr( randn(n,n) );
    [G,~] = qr( randn(m,n) );
    c = zeros(n,1);

    K   = CovarianceMatrix(eye(), 'C');
    C   = CovarianceMatrix(eye(m), 'C');

    kalman = UltimateKalman();


    j = 1;

    tic;

    for i=1:k
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

t = t(10:end,:);

r = (t(2:end,:) - t(1:end-1,:))/100;
r = r*1000; % milliseconds

medians = median(r)
size(r)
sizes

close all
figure
axis square
set(gca,'Box','on');;
hold on;
%plot(r,'k-','LineWidth',1);
plot(100*(0:size(r,1)-1),r);
xlim([0 100*(size(r,1)-1)]);
ylim([0 1.33*max(medians)])
legend(num2str(sizes'),'Location','SouthEast');
xlabel('step');
ylabel('time (ms)')
hold off;
%exportgraphics(gca,'../outputs/stress_6_96.pdf');



