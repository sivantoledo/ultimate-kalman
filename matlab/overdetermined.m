clear all; 

obsDim = 1; % number of rows in each observation, from 1 to 6

alpha = 2*pi/16;

F = [ cos(alpha) -sin(alpha) 
      sin(alpha)  cos(alpha) ];

G = [ 1 0 
      0 1 
      1 1
      2 1
      1 2
      3 1 ];

G = G(1:obsDim, :);

evolutionStd   = 1e-3;
observationStd = 1e-1;

Ce = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(2), 'C');
Co = CovarianceMatrix(observationStd*observationStd*eye(obsDim), 'C');

k = 16;

states = NaN * zeros(2,k);
obs    = NaN * zeros(obsDim,k);
states(:,1) = [ 1 ; 0 ];

for i=2:k
    states(:,i) = F*states(:,i-1) + evolutionStd * randn(2,1);
end

for i=1:k
    obs(:,i) = G*states(:,i) + observationStd * randn(obsDim,1);
end

kalman = UltimateKalman();
filtered = NaN * zeros(2,k);
for i=1:k
    kalman.advance(2,NaN);
    if (i>1); 
        kalman.evolve([],F,zeros(2,1),Ce);
    end
    kalman.observe(G,obs(:,i),Co);

    filtered(:,i) = kalman.filtered();
    % uncomment to test dropping (but then you can't smooth the entire
    % track.
    %kalman.drop();
end

kalman.smooth();

smoothed = NaN * zeros(2,k);
for i=1:k 
    smoothed(:,i) = kalman.smoothed(i-1); % step numbers are zero based
end

close all; 
figure;

axis square
set(gca,'Box','on');;
hold on;

plot(states(1,:),states(2,:),'k-','LineWidth',1);
if obsDim==2; plot(obs(1,:),obs(2,:),'r.'); end
plot(filtered(1,:),filtered(2,:),'b-','LineWidth',1);
plot(smoothed(1,:),smoothed(2,:),'m-','LineWidth',1);

xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
hold off;
exportgraphics(gca,'Underdetermined.pdf');





