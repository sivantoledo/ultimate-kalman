function projectile(kalman_factory,seed,exp)
% PROJECTILE a 2-dimensional projectile example for UltimateKalman
%
%    PROJECTILE(factory,seed,exp) 
%      factory: a handle to a function that returns UltimateKalman objects
%      seed:    random-number generator seed (default is 1)
%      exp:     whether to export the graphs to PDF files (default is false)
%
% copyright 2022-2024 Sivan Toledo

if nargin < 2
    seed = 1;
end

if nargin < 3
    exp = false;
end

rng(seed);

n = 4;
k = 1201;

delta_t = 0.1;
b = 1e-4;

F = [ 1  0  delta_t        0
      0  1        0  delta_t
      0  0    1 - b        0
      0  0        0    1 - b ];
G = [ 1 0 0 0
      0 1 0 0 ];

c = [ 0
      0 
      0 
      -9.8 * delta_t ];

l = size(F,1);
m = size(G,1);

evolutionStd   = sqrt(0.1);
observationStd = sqrt(5000);

K   = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(l), 'C');
C   = CovarianceMatrix(observationStd*observationStd*eye(m), 'C');

states = NaN * zeros(n,k);
obs    = NaN * zeros(m,k);

% simulate
states(:,1) = [ 0 
                0 
                300
                600 ]; % initial state
for i=2:k
   states(:,i) = F*states(:,i-1) + c + evolutionStd * randn(l,1);
end

% generate observations
for i=400:600
    obs(:,i) = G*states(:,i) + observationStd * randn(m,1);
end

% we use observations 400 to 600 only.
% 
% First, generate an initial expectation for a traditional
% Kalman filter.

initialExpectation = [ obs(1,400) 
                       obs(2,400)
                       (obs(1,410) - obs(1,400)) / (10*delta_t) 
                       (obs(2,410) - obs(2,400)) / (10*delta_t) ];

kalman = kalman_factory();
filtered = NaN * zeros(n,k);
smoothed = NaN * zeros(n,k);
est      = NaN * zeros(n,k);
std      = NaN * zeros(n,k);
filtstd  = NaN * zeros(n,k);
smthstd  = NaN * zeros(n,k);

kalman.evolve(4);
kalman.observe(eye(4),initialExpectation,CovarianceMatrix(1e5*eye(4),'C'));

% evolve and observe
for i=401:600
    kalman.evolve(4,[],F,c,K);
    kalman.observe(G,obs(:,i),C);
    filtered(:,i) = kalman.estimate();
end

% continue to predict

% evolve and observe
for i=601:k
    kalman.evolve(4,[],F,c,K);
    kalman.observe();
    filtered(:,i) = kalman.estimate();
end

% smooth to get point of origin
% now we start the filter from 0, but 
% we do not provide the initial expectation
kalman = kalman_factory();
% evolve and observe
for i=1:k
    kalman.evolve(4,[],F,c,K);
    if i >= 400 && i <= 600
        kalman.observe(G,obs(:,i),C);
    else
        kalman.observe();
    end
end

kalman.smooth();

% collect smoothed estimates
for i=1:k
    smoothed(:,i) = kalman.estimate(i-1);
end


%close all
figure
axis square
set(gca,'Box','on');;
hold on;
plot(states(1,:),states(2,:),'k-','LineWidth',1);
plot(smoothed(1,:),smoothed(2,:),'m-','LineWidth',1);
plot(filtered(1,:),filtered(2,:),'b-','LineWidth',1);
plot(obs(1,:),obs(2,:),'r.');
%plot(full(:,1),full(:,2),'m-','LineWidth',1);
hold off;
if exp; exportgraphics(gca,'../outputs/projectile.pdf'); end;

figure
axis square
set(gca,'Box','on');;
hold on;
plot(states(1,:),states(2,:),'k-','LineWidth',1);
plot(obs(1,:),obs(2,:),'r.');
plot(filtered(1,:),filtered(2,:),'b-','LineWidth',1);
%plot(full(:,1),full(:,2),'m-','LineWidth',1);
xlim([ min(obs(1,:)) max(obs(1,:)) ]);
ylim([ min(obs(2,:)) max(obs(2,:)) ]);
hold off;
if exp; exportgraphics(gca,'../outputs/projectile_zoom.pdf'); end;



