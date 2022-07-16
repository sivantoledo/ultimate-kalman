function add_remove(seed)

if nargin<1
    seed = 1;
end

rng(seed);

n = 2;
k = 10;

state = [ 1
          2 ];

% selection matrices

S_1    = [ 1 0 ];
S_2    = [ 0 1 ];
S_both = [ 1 0
           0 1 ];

F1 = eye(1);
F2 = eye(2);
G1 = eye(1);
G2 = eye(2);
c = [ 0 ];

evolutionStd   = 0.1;
observationStd = 0.1;

K1   = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(1), 'C');
C1   = CovarianceMatrix(observationStd*observationStd*eye(1), 'C');
K2   = CovarianceMatrix(evolutionStd  *evolutionStd  *eye(2), 'C');
C2   = CovarianceMatrix(observationStd*observationStd*eye(2), 'C');

c1 = zeros(1,1);
c2 = zeros(2,1);

kalman = UltimateKalman();
obs      = NaN * zeros(n,k);

% generate observations
for i=1:k
    obs(:,i) = state + observationStd * randn(n,1);
end

% step 0
kalman.evolve(1);
kalman.observe(G1,S_1*obs(:,1),C1);
e = kalman.estimate()

% step 1
kalman.evolve(1,[],F1,c1,K1);
kalman.observe(G1,S_1*obs(:,2),C1);
e = kalman.estimate()

% step 2: add the second parameter using H
kalman.evolve(2,[ 1 ; 0],F1,c1,K1);
kalman.observe(G2,S_both*obs(:,3),C2);
e = kalman.estimate()

% step 3
kalman.evolve(2,[],F2,c2,K2);
kalman.observe(G2,S_both*obs(:,4),C2);
e = kalman.estimate()

% step 4: drop back to one parameter, the second; select using F
kalman.evolve(1,[],[0 1],c1,K1);
kalman.observe(G1,S_2*obs(:,5),C1);
e = kalman.estimate()

% step 5
kalman.evolve(1,[],F1,c1,K1);
kalman.observe(G1,S_2*obs(:,1),C1);
e = kalman.estimate()

