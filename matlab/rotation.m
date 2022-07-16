function rotation(seed,obs_dim)

if nargin<2
    obs_dim = 2; % number of rows in each observation, from 1 to 6
end

if nargin>=1
    rng(seed);
end

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

K = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(2),      'C');
C = CovarianceMatrix(observationStd*observationStd*eye(obs_dim), 'C');

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

%[states' obs']

kalman = UltimateKalman();
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
    kalman.observe(G);
    predicted(:,i) = kalman.estimate(i-1); % zero based step numbers
end

% roll back and add observation of step 1
kalman.rollback(1);
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

kalman.smooth();

smoothed = NaN * zeros(2,k);
for i=1:k 
    smoothed(:,i) = kalman.estimate(i-1); % zero based step numbers
end

close all; 
figure;

axis square
set(gca,'Box','on');
hold on;

plot(states(1,:),states(2,:),'k-','LineWidth',1);
if obs_dim==2; plot(obs(1,:),obs(2,:),'r.'); end
plot(filtered(1,:),filtered(2,:),'b-','LineWidth',1);
plot(smoothed(1,:),smoothed(2,:),'m-','LineWidth',1);
plot(predicted(1,:),predicted(2,:),'c-','LineWidth',1);

xlim([-1.25 1.25]);
ylim([-1.25 1.25]);
hold off;
%exportgraphics(gca,'Underdetermined.pdf');
