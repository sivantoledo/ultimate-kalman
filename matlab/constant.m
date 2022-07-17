function constant(seed, k, smooth, sloping, exception)
% constant(seed, k, smooth, sloping, exception) a test function for UltimateKalman
%
% copyright 2022 Sivan Toledo

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

kalman = UltimateKalman();
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

close all
figure
axis square
set(gca,'Box','on');;
hold on;
plot(t,states,'k-','LineWidth',1);
if (smooth)
  plot(t,est,'m-','LineWidth',1);
else
  plot(t,est,'b-','LineWidth',1);
end

xs = [t, fliplr(t)];
ys = [est+3*std, fliplr(est-3*std)];
fill(xs, ys,'k','FaceAlpha',0.2,'LineStyle','none');

%obs
%min(obs(:,1))
%max(obs(:,1));
%states

plot(t,obs,'r.');

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


