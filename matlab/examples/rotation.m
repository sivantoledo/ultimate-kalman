function rotation(kalman_factory,seed,obs_dim)
% ROTATION a simple test function for UltimateKalman, tracking a rotating 
%          point in the plane
%
%    ROTATION(factory, seed, obs_dim) a test function for UltimateKalman
%      factory: a handle to a function that returns UltimateKalman objects
%      seed:    random-number generator seed
%      obs_dim: the number of observations per time step, between 1 and 6
%
% copyright 2022-2024 Sivan Toledo

rng(seed);

alpha = 2*pi/16;

F = [ cos(alpha) -sin(alpha) 
      sin(alpha)  cos(alpha) ];

G = [ 1 0 
      0 1 
      1 1
      2 1
      1 2
      3 1 ];

G = G(1:obs_dim, :);

evolutionStd   = 1e-3;
observationStd = 1e-1;

K = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(2),       'C');
C = CovarianceMatrix(observationStd*observationStd*eye(obs_dim), 'C');
K = CovarianceMatrix(evolutionStd^-1  *ones(2,1),'w');
C = CovarianceMatrix(observationStd^-1*ones(obs_dim,1), 'w');

k = 16;

states = NaN * zeros(2,k);
obs    = NaN * zeros(obs_dim,k);
states(:,1) = [ 1 ; 0 ];

for i=2:k
    states(:,i) = F*states(:,i-1) + evolutionStd * randn(2,1);
end

for i=1:k
    obs(:,i) = G*states(:,i) + observationStd * randn(obs_dim,1);
end

states
obs

%[states' obs']

kalman = kalman_factory();
filtered  = NaN * zeros(2,k);
predicted = NaN * zeros(2,k);

if true
% do the first step
kalman.evolve(2,[],F,zeros(2,1),K);
kalman.observe(G,obs(:,1),C);
predicted(:,1) = kalman.estimate(0); 
% predict
for i=2:k
    kalman.evolve(2,[],F,zeros(2,1),K);
    kalman.observe();
    predicted(:,i) = kalman.estimate(i-1); % zero based step numbers
end
predicted
if true
% roll back and add observation of step 1
kalman.earliest
kalman.latest
kalman.rollback(1);
disp('rollback')
kalman.earliest
kalman.latest
%kalman.evolve(2,[],F,zeros(2,1),K);
kalman.observe(G,obs(:,2),C);
end

filtered(:,1) = kalman.estimate(0);
filtered(:,2) = kalman.estimate(1);
% evolve and observe steps 2 and up
for i=3:k
    kalman.evolve(2,[],F,zeros(2,1),K);

    kalman.observe(G,obs(:,i),C);

    filtered(:,i) = kalman.estimate(i-1); % zero based step numbers
    % uncomment to test dropping (but then you can't smooth the entire
    % track.
    %kalman.forget();
end

filtered
end

kalman.smooth();

E = [];

smoothed = NaN * zeros(2,k);
for i=1:k 
    [est,cov] = kalman.estimate(i-1); % zero based step numbers
    smoothed(:,i) =  est;
    if i==1 
        E = cov;
    end
end

smoothed

disp('covariance matrix of smoothed state 0:');
Eexplicit = E.explicit()
disp('W such that W^T*W = cov^-1:');
W = E.weigh(eye(2))

fprintf('std deviation of first state coordinate %.2e (of first observation %.2e)\n',...
        sqrt(Eexplicit(1,1)),observationStd);

gray   = [7.3333e-01   7.3333e-01   7.3333e-01];
red    = [9.3333e-01   4.0000e-01   4.6667e-01];
cyan   = [4.0000e-01   8.0000e-01   9.3333e-01];
green  = [1.3333e-01   5.3333e-01   2.0000e-01];
blue   = [2.6667e-01   4.6667e-01   6.6667e-01];
yellow = [8.0000e-01   7.3333e-01   2.6667e-01];

linewidth  = 1.25;
markersize = 8;

%close all; 
figure;

axis square
set(gca,'Box','on');
hold on;

plot(states(1,:),states(2,:),'Color',yellow,'LineWidth',linewidth,'Marker','none');
if obs_dim==2; 
    plot(obs(1,:),obs(2,:),'Marker','.','MarkerSize',markersize,'LineStyle','none','Color',red); 
end
plot(filtered(1,:),filtered(2,:),'Color',blue,'LineWidth',linewidth,'Marker','none');
plot(smoothed(1,:),smoothed(2,:),'Color',cyan,'LineWidth',linewidth,'Marker','none');
plot(predicted(1,:),predicted(2,:),'Color',green,'LineWidth',linewidth,'Marker','none');

xlim([-1.25 1.25]);
ylim([-1.25 1.25]);
hold off;
%exportgraphics(gca,'Underdetermined.pdf');
