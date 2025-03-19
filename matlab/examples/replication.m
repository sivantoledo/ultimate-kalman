function replication(kalman_factory,exp,perf)
% REPLICATION testing of an UltimateKalman implementation and optionally
%             export graphs and test performance
%
%    REPLICATION(kalman_factory,exp,perf) 
%      kalman_factory: a factory function returned by kalmanFactory
%      exp:            boolean, whether to export graphs (default is false)
%      perf:           an optional cell array of Kalman factories whose
%                      performance the function should measure and plot
%
% copyright 2022-2024 Sivan Toledo

if nargin<1
    kalman_factory = kalmanFactory('KalmanUltimate'); 
    perf = {kalmanFactory('KalmanUltimate'),kalmanFactory('KalmanJava'),kalmanFactory('KalmanNative')};
end

if nargin<2
    exp = false;
end

if nargin<3 && nargin>=1
    perf = {};
end

if exp
    mkdir '../../outputs';
end

if length(perf)>0
performance(perf, 1 , [6],1e5,1000);
if exp; exportgraphics(gca,'../../outputs/perftest_imps_6.pdf');  end;

performance(perf, 1 , [48],1e5,1000);
if exp; exportgraphics(gca,'../../outputs/perftest_imps_48.pdf'); end;

%performance({'C'}, 1 , [6 12 24 48 96],1e5,1000);
%if exp; exportgraphics(gca,'../../outputs/perftest_C_6_96.pdf');  end;
end

clock_offsets(kalman_factory, 6)
if exp; exportgraphics(gca,'../../outputs/clock_offsets.pdf'); end;

rotation(kalman_factory, 5,2)
if exp; exportgraphics(gca,'../../outputs/rotation2.pdf'); end;
rotation(kalman_factory, 5,1)
if exp; exportgraphics(gca,'../../outputs/rotation1.pdf'); end;

% constant(kalman_factory, 1,101,true,false,[NaN,1])
% if exp; exportgraphics(gca,'../../outputs/constant_smoothed.pdf'); end;
% constant(kalman_factory, 1,101,false,false,[NaN,1])
% if exp; exportgraphics(gca,'../../outputs/constant_filtered.pdf'); end;

constant(kalman_factory, 1,101,true,true,[NaN,1])
if exp; exportgraphics(gca,'../../outputs/sloping_smoothed.pdf'); end;
constant(kalman_factory, 1,101,false,true,[NaN,1])
if exp; exportgraphics(gca,'../../outputs/sloping_filtered.pdf'); end;

constant(kalman_factory, 1,101,true,true,[50,0.25])
if exp; exportgraphics(gca,'../../outputs/sloping_smoothed_exception.pdf'); end;
constant(kalman_factory, 1,101,false,true,[50,0.25])
if exp; exportgraphics(gca,'../../outputs/sloping_filtered_exception.pdf'); end;

% projectile
% this function saves the plots to pdf files, so no
% need to call exportgraphics here.
projectile(kalman_factory,1,exp);

% The next script should produce two steps with an estimate near 1,
% than 2 more with one estimate near 1 and the other near 2, and then
% the first parameter is removed, so values near 1 should be returned
% twice more.
add_remove(kalman_factory,1);

end