function clock_offsets(kalman_factory,seed)
% CLOCK_OFFSETS a test function for UltimateKalman, demonstrating how to 
%               track clock offsets in a distributed system. The example
%               also shows how to model parameters that affect only a
%               single time step.
%
%    CLOCK_OFFSETS(factory, seed) a test function for UltimateKalman
%      factory: a handle to a function that returns UltimateKalman objects
%      seed:    random-number generator seed
%
% copyright 2022-2024 Sivan Toledo
%
% copyright 2022 Sivan Toledo

if nargin < 1
    seed = 4;
end

rng(seed);

k = 251;
clockCount = 3;
evolutionStd   =  10e-9;
observationStd = 100e-9;
initialOffsetStd = 1e-6;
delayStd         = 10e-6;

F  = [eye(clockCount) zeros(clockCount,1) ];
H  = [eye(clockCount) zeros(clockCount,1) ];
G  = [eye(clockCount)  ones(clockCount,1) ];
G0 = [eye(clockCount)  ones(clockCount,1) 
                   1  zeros(1,clockCount) ];
c  = zeros(clockCount,1);

n = clockCount + 1;
l = size(F,1);
m = size(G,1);

K   = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(l)  , 'C');
C   = CovarianceMatrix(observationStd*observationStd*eye(m)  , 'C');
C0  = CovarianceMatrix(observationStd*observationStd*eye(m+1), 'C');

%%%%%% Simulation %%%%%%

t    = NaN * zeros(k, clockCount); % arrival times
tau  = (0:k-1)';                   % departure times; exact but we do not use this information
d    = delayStd         *   randn(1,clockCount);
f    = initialOffsetStd * [ randn(1, clockCount) ; zeros(k-1,clockCount) ]; % the offsets will evolve; zeros are for memory allocation
r    = evolutionStd     *   randn(1,clockCount);    % rate errors 

for i=1:k
    for j=1:clockCount
        if (i>1) % initial conditions defined above
            f(i,j) = f(i-1,j) + evolutionStd * randn;
            %o(i,j) = o(i-1,j) + r(j);
        end
        t(i,j) = tau(i) + d(j) + f(i,j) + observationStd * randn;
    end
end

kalman = kalman_factory();
filtered = NaN * zeros(n,k);

kalman.evolve(n);
kalman.observe(G0, [(t(1,:)-d) 0]', C0);
filtered(:,1) = kalman.estimate();

% evolve and observe
for i=2:k
    kalman.evolve(n,H,F,c,K);
    kalman.observe(G,(t(i,:)-d)',C);
    filtered(:,i) = kalman.estimate();
end

% 

relative = f - f(:,1);

%close all;
figure
axis square
set(gca,'Box','on');;
hold on;
plot(tau',relative);
set(gca,'ColorOrderIndex',1);
plot(tau,filtered(1:clockCount,:)-filtered(1,:),'--');
xlabel('time (s)');
ylabel('relative offset (s)');
hold off;
%exportgraphics(gca,'../outputs/projectile_zoom.pdf');



