function kalmanTest()

stateThreshold = 1e-8;
covThreshold   = 1e-4;

rng(1)

% k = KalmanUltimate();
% k = KalmanNative(struct('algorithm','Ultimate'));
% n = 3;
% k.evolve(n);
% Cf = [ 4 -1 -1
%        -1 3 -1
%        -1 -1 3]
% L = chol(Cf,'lower')
% W = inv(L)
% CC = CovarianceMatrix(Cf,'C');
% CW = CovarianceMatrix(W,'W');
% k.observe(eye(n),[1 2 3]',CC);
% [e,c] = k.estimate();
% e
% c
% W=c.W
% inv(W'*W)
% c.C
% return


NONE          = 0;
UNOBSERVABALE = 1;     
RECT_H        = 2;    
LIMITATION_PLURAL = 4; % a restriction: only works on k>1    

msparse = { 'msparse', kalmanFactory('KalmanSparse'),                                                bitor( UNOBSERVABALE, RECT_H) };
mconv   = { 'mconv',   kalmanFactory('KalmanConventional'),                                          0 };
nconv   = { 'nconv',   kalmanFactory('KalmanNative',struct('algorithm','Conventional')),             0 };
pconv   = { 'pconv',   kalmanFactory('KalmanPython',struct('algorithm','KalmanConventional')),       0 };
multps  = { 'multps',  kalmanFactory('KalmanUltimate',struct('estimateCovariance','PaigeSaunders')), bitor( UNOBSERVABALE, RECT_H) };
pultps  = { 'pultps',  kalmanFactory('KalmanPython',struct('algorithm','KalmanUltimate')),           bitor( UNOBSERVABALE, RECT_H) };
multsi  = { 'multsi',  kalmanFactory('KalmanUltimate',struct('estimateCovariance','SelInv')),        bitor( UNOBSERVABALE, RECT_H) };
nultps  = { 'nultps',  kalmanFactory('KalmanNative',struct('algorithm','Ultimate')),                 bitor( UNOBSERVABALE, RECT_H) };
jultps  = { 'jultps',  kalmanFactory('KalmanJava'),                                                  bitor( UNOBSERVABALE, RECT_H) };
moddevn = { 'moddevn', kalmanFactory('KalmanOddevenSmoother',struct()),          bitor(bitor( UNOBSERVABALE, RECT_H), LIMITATION_PLURAL) };
noddevn = { 'noddevn', kalmanFactory('KalmanNative',struct('algorithm','Oddeven','estimateCovariance',false)),                  bitor(bitor( UNOBSERVABALE, RECT_H), LIMITATION_PLURAL) };
massoc  = { 'massoc',  kalmanFactory('KalmanAssociativeSmoother'),                                   LIMITATION_PLURAL  };

variants = { msparse
             moddevn
             massoc
             multps
             nultps 
             mconv
             nconv
           };

variants = { msparse
             multps
             multsi
             pultps
             nultps
             mconv
             pconv
             nconv
             moddevn
             noddevn
             massoc
           };

variants = { multps
             multsi
             nultps
             moddevn
             noddevn
           };

dimensions = { uniform(1,1,1,1)
          uniform(1,1,1,7)
          uniform(2,2,2,1)
          uniform(2,2,2,7)
          uniform(3,3,3,1)
          uniform(3,3,3,4)
          uniform(3,3,3,5)
          uniform(3,3,1,3) 
          uniform(4,4,2,2) 
          uniform(2,1,2,9)
          uniform(3,2,3,9)
          uniform(3,1,3,9)
          [3 0 2 nan nan ; 3 1 2 nan nan ; 3 1 3 nan nan ; 3 1 2 nan nan ] % note the state with 3 observation in the middle; it's necessary
          [1 1 0 nan nan ; 1 1 2 nan nan   ] % increasing state dimension
          [1 1 0 nan nan ; 1 1 1 nan nan ; 2 1 2 nan nan ; 3 2 3 nan nan ; 3 3 0 nan nan ; 3 2 3 nan nan ; 4 2 3 nan nan] % increasing state dimension
        }

%dimensions = { uniform(1,1,1,1) }
%dimensions = {[3 0 2 nan nan ; 3 1 2 nan nan ; 3 1 3 nan nan ; 3 1 2 nan nan ]}; % high covariance errors, ~ 1e-7

%          uniform(2,1,1,9)

[v,FS,FC,SS,SC] = testDimensions(dimensions);
v
FS
FC
SS
SC

%return

[v,FS,FC,SS,SC,covTypes] = testCovarianceTypes(uniform(9,9,9,17),{ 'W', 'C', 'w' });
covTypes
v
FS
FC
SS
SC

function [v,FS,FC,SS,SC] = testDimensions(tests)

n = length(variants)-1;
m = size(tests,1);

FS=NaN*ones(n,m);
FC=NaN*ones(n,m);
SS=NaN*ones(n,m);
SC=NaN*ones(n,m);

v={};

for i=1:n
    for j=1:m
        %[i j]
        kalman1 = variants{1}{2}();
        kalman2 = variants{i+1}{2}();
        capabilities1 = variants{1}{3};
        capabilities2 = variants{i+1}{3};
        %class(kalman1)
        %class(kalman2)
        test = tests{j,1};
        if rect_h(test)     && ~bitand(capabilities1,RECT_H       ); continue; end
        if observable(test) && ~bitand(capabilities1,UNOBSERVABALE); continue; end
        if rect_h(test)     && ~bitand(capabilities2,RECT_H       ); continue; end
        if observable(test) && ~bitand(capabilities2,UNOBSERVABALE); continue; end
        if ~plural(test)     && bitand(capabilities2,LIMITATION_PLURAL); continue; end
        test
        v{i,1}=variants{i+1}{1}
        %plural(test)
        [fse,fce,sse,sce] = compare(kalman1, kalman2, test, 'W', 'W');
        if isfinite(fse) && fse  > stateThreshold; FS(i,j) = 0; end 
        if isfinite(fce) && fce  > covThreshold; FC(i,j) = 0; end 
        if isfinite(sse) && sse  > stateThreshold; SS(i,j) = 0; end 
        if isfinite(sce) && sce  > covThreshold; SC(i,j) = 0; end 
        if isfinite(fse) && fse <= stateThreshold; FS(i,j) = 1; end 
        if isfinite(fce) && fce <= covThreshold; FC(i,j) = 1; end 
        if isfinite(sse) && sse <= stateThreshold; SS(i,j) = 1; end 
        if isfinite(sce) && sce <= covThreshold; SC(i,j) = 1; end 
        if ~isfinite(fse); FS(i,j) = fse; end 
        if ~isfinite(fce); FC(i,j) = fce; end 
        if ~isfinite(sse); SS(i,j) = sse; end 
        if ~isfinite(sce); SC(i,j) = sce; end 
    end
    
end

end

function [v,FS,FC,SS,SC,covTypes] = testCovarianceTypes(test,covTypes)
%%% test covariance matrix types

%test = uniform(1,1,1,17);

n = length(variants);
m = length(covTypes);

FS=NaN*ones(n,m);
FC=NaN*ones(n,m);
SS=NaN*ones(n,m);
SC=NaN*ones(n,m);

v={};

for i=1:n
    for j=1:m
        %[i j]
        %kalman1 = multps{2}();
        kalman1 = msparse{2}();
        kalman2 = variants{i}{2}();
        v{i,1}=variants{i}{1}
        [fse,fce,sse,sce] = compare(kalman1, kalman2, test, 'W', covTypes{j});
        if isfinite(fse) && fse  > stateThreshold; FS(i,j) = 0; end 
        if isfinite(fce) && fce  > covThreshold; FC(i,j) = 0; end 
        if isfinite(sse) && sse  > stateThreshold; SS(i,j) = 0; end 
        if isfinite(sce) && sce  > covThreshold; SC(i,j) = 0; end 
        if isfinite(fse) && fse <= stateThreshold; FS(i,j) = 1; end 
        if isfinite(fce) && fce <= covThreshold; FC(i,j) = 1; end 
        if isfinite(sse) && sse <= stateThreshold; SS(i,j) = 1; end 
        if isfinite(sce) && sce <= covThreshold; SC(i,j) = 1; end 
        if ~isfinite(fse); FS(i,j) = fse; end 
        if ~isfinite(fce); FC(i,j) = fce; end 
        if ~isfinite(sse); SS(i,j) = sse; end 
        if ~isfinite(sce); SC(i,j) = sce; end 
    end
    
end

end

function prop = rect_h(t)
    if sum(t(:,1) > t(:,2)) > 0 % n > l for some row
        prop = true;
    else
        prop = false;
    end
end

function prop = observable(t)
    if sum(t(1,1) > t(1,3)) > 0 % n > m in first row
        prop = true;
    else
        prop = false;
    end
end

function prop = plural(t)
    if size(t,1) > 1 % more than 1 step
        prop = true;
    else
        prop = false;
    end
end


function t = uniform(n,l,m,k)
t = ones(k,1) * [ n l m nan nan];
end

function t = predict(n,l,m,k)
t = ones(k,1) * [ n l 0 nan nan];
t(1,:) = [ n 0 m nan nan ];
end

end