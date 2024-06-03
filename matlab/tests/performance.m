function t = performance(kalmans, seed, sizes, count, decimation)
% PERFORMANCE a function to test the performance of UltimateKalman
%             implementations
%
%    PERFORMANCE(factories, seed, sizes, count, decimation)
%      factories:  a cell array of functions that return UltimateKalman
%                  objects
%      seed:       random-number generator seed
%      sizes:      state-vector dimensions to test the implemenations on 
%      count:      how many time steps to simulate and filter
%      decimation: number of time steps between timestamps
%
% copyright 2022-2024 Sivan Toledo

if nargin < 1
    seed = 1;
end

rng(seed);

q = 1;

for k = 1:length(kalmans)
    implementation = kalmans{k}

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

    kalman = factory(implementation);

        tstart=tic;
        t(:,q) = kalman.perftest([],F,c,K,G,o,C,count,decimation);
        tend=toc(tstart);
        ave_ms = 1000*tend/count
        q = q+1;
    end
end

%t = t(10:end,:);

%r = (t(2:end,:) - t(1:end-1,:))/100;
r = 1000*t; % r*1000; % milliseconds

%r

medians = median(r)
size(r)
sizes

%close all
figure
axis square
set(gca,'Box','on');;
hold on;
%plot(r,'k-','LineWidth',1);
plot(decimation*(0:size(r,1)-1),r,'LineWidth',1);
xlim([0 decimation*(size(r,1)-1)]);
ylim([0 1.33*max(medians)])
if length(sizes) > 1
  legend(num2str(sizes'),'Location','SouthEast');
  title(kalmans{1})
else
  legend(kalmans,'Location','SouthEast');
  title(sprintf('n = %d',sizes(1)));
end
xlabel('step');
ylabel('time (ms)')
hold off;
%exportgraphics(gca,'../outputs/stress_6_96.pdf');

function kalman = factory(implementation)
    switch implementation
        case 'C'
            kalman = UltimateKalmanNative;
        case 'Java'
            kalman = UltimateKalmanJava;
        case 'MATLAB'
            kalman = UltimateKalman;
        otherwise
            error("implementation should be 'M', 'N', or 'J' (for MATLAB, native, or Java)");
    end
end

end



