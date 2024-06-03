function replication(implementation,exp,perf)
% REPLICATION replication of testing of UltimateKalman
%
%    REPLICATION(implementation,exp,perf) performs tests
%      on a given implementation ('MATLAB', 'C', or 'Java')
%      If exp==true, graphs are exported to PDF files.
%      if perf==true, the function also runs performance tests.
%      The defaults are 'MATLAB', no export, and no performance tests.
%
% copyright 2022-2024 Sivan Toledo

if nargin<1
    implementation = 'MATLAB'; % matlab
end

if nargin<2
    exp = false;
end

if nargin<3
    perf = false;
end

addpath '..';

if (~isempty(ver('MATLAB')))
    disp('running under MATLAB, adding native and java support')
    addpath '..\..\native'
    javaaddpath('../../java/ultimatekalman.jar');
    javaaddpath('../../java/commons-math3-3.6.1.jar');
end

if (~isempty(ver('Octave')))
    disp('running under Octave, adding native and java support')
    addpath '../../..\native'
    javaaddpath('../../ultimatekalman.jar');
    javaaddpath('../../commons-math3-3.6.1.jar');
end

if perf
performance({'MATLAB', 'Java', 'C'}, 1 , [6],1e5,1000);
if exp; exportgraphics(gca,'../../outputs/perftest_imps_6.pdf');  end;

performance({'MATLAB', 'Java', 'C'}, 1 , [48],1e5,1000);
if exp; exportgraphics(gca,'../../outputs/perftest_imps_48.pdf'); end;

%performance({'C'}, 1 , [6 12 24 48 96],1e5,1000);
%if exp; exportgraphics(gca,'../../outputs/perftest_C_6_96.pdf');  end;
end

clock_offsets(@factory, 6)
if exp; exportgraphics(gca,'..\..\outputs\clock_offsets.pdf'); end;

rotation(@factory, 5,2)
if exp; exportgraphics(gca,'..\..\outputs\rotation2.pdf'); end;
rotation(@factory, 5,1)
if exp; exportgraphics(gca,'..\..\outputs\rotation1.pdf'); end;

constant(@factory, 1,101,true,false,[NaN,1])
if exp; exportgraphics(gca,'..\..\outputs\constant_smoothed.pdf'); end;
constant(@factory, 1,101,false,false,[NaN,1])
if exp; exportgraphics(gca,'..\..\outputs\constant_filtered.pdf'); end;

constant(@factory, 1,101,true,true,[50,0.25])
if exp; exportgraphics(gca,'..\..\outputs\sloping_smoothed_exception.pdf'); end;
constant(@factory, 1,101,false,true,[50,0.25])
if exp; exportgraphics(gca,'..\..\outputs\sloping_filtered_exception.pdf'); end;

% projectile
% this function saves the plots to pdf files, so no
% need to call exportgraphics here.
projectile(@factory,1,exp);

% The next script should produce two steps with an estimate near 1,
% than 2 more with one estimate near 1 and the other near 2, and then
% the first parameter is removed, so values near 1 should be returned
% twice more.
add_remove(@factory,1);

function kalman = factory()
    switch implementation
        case 'C'
            kalman = UltimateKalmanNative;
        case 'Java'
            kalman = UltimateKalmanJava;
        case 'MATLAB'
            kalman = UltimateKalman;
        otherwise
            error("implementation should be 'MATLAB', 'C', or 'Java'");
    end
end

end