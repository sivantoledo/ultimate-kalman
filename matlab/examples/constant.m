function constant(kalman_factory, seed, k, smooth, sloping, exception)
% CONSTANT a simple test function for UltimateKalman; simulates a scalar
%          target whose dynamics are either staying constant or increasing
%          at a fixed rate
%
%    CONSTANT(factory, seed, k, smooth, sloping, exception) a test function for UltimateKalman
%      factory: a handle to a function that returns UltimateKalman objects
%      seed:    random-number generator seed
%      k:       number of steps to simulate and filter
%      smooth:  boolean, whether to smooth or just filter
%      sloping: boolean, whether the evolution rule should cause the target
%               to rise (slope up) or stay constant; evolution is with
%               noise either way
%      exception: vector with step-number for exceptional step and standard
%               deviation for this step (first and second positions)
%
% copyright 2022-2024 Sivan Toledo

rng(seed);

n = 1;

F = [ 1 ];
G = [ 1 ];
c = [ 0 ];

l = size(F,1);
m = size(G,1);

evolutionStd   =  1;
observationStd = 10;
exceptionStd   = exception(2);

K   = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(l), 'C');
C   = CovarianceMatrix(observationStd*observationStd*eye(m), 'C');
Cex = CovarianceMatrix(exceptionStd  *exceptionStd  *eye(m), 'C');


K   = CovarianceMatrix(evolutionStd^-1  *ones(1,1), 'w');
C   = CovarianceMatrix(observationStd^-1*ones(m,1), 'w');

states = NaN * zeros(1,k);
obs    = NaN * zeros(1,k);

% simulate
states(:,1) = [ 0 ]; % initial state
for i=2:k
    if (sloping)
      states(:,i) = F*states(:,i-1) + 0.2;
    else
      states(:,i) = F*states(:,i-1) + evolutionStd * randn(l,1);
    end
end

% generate observations
for i=1:k
    if  i == exception(1)
        obs(:,i) = G*states(:,i) + exceptionStd   * randn(m,1);
    else
        obs(:,i) = G*states(:,i) + observationStd * randn(m,1);
    end
end

kalman = kalman_factory();
filtered = NaN * zeros(n,k);
smoothed = NaN * zeros(n,k);
est      = NaN * zeros(n,k);
std      = NaN * zeros(n,k);
filtstd  = NaN * zeros(n,k);
smthstd  = NaN * zeros(n,k);

% evolve and observe
for i=1:k
    kalman.evolve(1,[],F,c,K);

    if  i == exception(1)
        kalman.observe(G,obs(:,i),Cex);
    else
        kalman.observe(G,obs(:,i),C);
    end
    
end

if smooth
    kalman.smooth();
end

for i=1:k
    %[ i kalman.earliest kalman.latest ]
    [est(:,i),cov] = kalman.estimate(i-1); % zero based step numbers
    std(:,i) = sqrt( cov.explicit() );
    kalman.forget(i-1);
end

t = 0:k-1;

% pallete = [ 2.6667e-01   4.6667e-01   6.6667e-01
%             4.0000e-01   8.0000e-01   9.3333e-01
%             1.3333e-01   5.3333e-01   2.0000e-01
%             8.0000e-01   7.3333e-01   2.6667e-01
%             9.3333e-01   4.0000e-01   4.6667e-01
%             7.3333e-01   7.3333e-01   7.3333e-01];

gray   = [7.3333e-01   7.3333e-01   7.3333e-01];
red    = [9.3333e-01   4.0000e-01   4.6667e-01];
cyan   = [4.0000e-01   8.0000e-01   9.3333e-01];
green  = [1.3333e-01   5.3333e-01   2.0000e-01];
blue   = [2.6667e-01   4.6667e-01   6.6667e-01];
yellow = [8.0000e-01   7.3333e-01   2.6667e-01];

linewidth  = 1.5;
markersize = 8;

% plot(1:100,1*ones(100,1),'LineWidth',2,'Color',gray);
% hold on
% plot(1:100,2*ones(100,1),'LineWidth',2,'Color',red);
% plot(1:100,3*ones(100,1),'LineWidth',2,'Color',cyan);
% plot(1:100,4*ones(100,1),'LineWidth',2,'Color',green);
% plot(1:100,5*ones(100,1),'LineWidth',2,'Color',blue);
% plot(1:100,6*ones(100,1),'LineWidth',2,'Color',yellow);
% ylim([0 7])
% hold off


%close all
figure
axis square
set(gca,'Box','on');
hold on;

xs = [t, fliplr(t)];
ys = [est+3*std, fliplr(est-3*std)];
%fill(xs, ys,'k','FaceAlpha',0.2,'LineStyle','none');
fill(xs, ys,gray,'LineStyle','none');


%plot(t,states,'k-','LineWidth',1);
h=plot(t,states,'Color',yellow,'LineWidth',linewidth,'Marker','none');
if (smooth)
  plot(t,est,'Color',cyan,'LineWidth',linewidth,'Marker','none');
else
  plot(t,est,'Color',blue,'LineWidth',linewidth,'Marker','none');
end


%obs
%min(obs(:,1))
%max(obs(:,1));
%states

%h=plot(t,obs,'r.');
%h
h=plot(t,obs,'Marker','.','MarkerSize',markersize,'LineStyle','none','Color',red);

xlim([ min(t) max(t) ]);
%ylim([ min(obs) max(obs) ]);

xlabel('step');

hold off;

%constant(1,101,false,true)
%exportgraphics(gca,'..\outputs\slope_filtered.pdf');
%constant(1,101,true,true)
%exportgraphics(gca,'..\outputs\slope_smoothed.pdf');
%constant(1,101,true,false)
%exportgraphics(gca,'..\outputs\constant_smoothed.pdf');
%constant(1,101,false,false)
%exportgraphics(gca,'..\outputs\constant_filtered.pdf');

%constant(1,101,true,true,[50,0.25])
%exportgraphics(gca,'..\outputs\sloping_smoothed_exception.pdf');
%constant(1,101,false,true,[50,0.25])
%exportgraphics(gca,'..\outputs\sloping_filtered_exception.pdf');


