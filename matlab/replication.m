function replication(implementation)
% test-replication script for UltimateKalman
%
% copyright 2022 Sivan Toledo

if nargin<1
    implementation = 'M'; % matlab
end

exp = false;

if (~isempty(ver('MATLAB')))
    disp('running under MATLAB, adding native and java support')
    exp = true
    addpath '..\native'
    javaaddpath('../ultimatekalman.jar');
    javaaddpath('../commons-math3-3.6.1.jar');
end

%performance({'MATLAB', 'Java', 'C'}, 1 , [6],1e5,1000);
%exportgraphics(gca,'../outputs/perftest_imps_6.pdf');

%performance({'MATLAB', 'Java', 'C'}, 1 , [48],1e5,1000);
%exportgraphics(gca,'../outputs/perftest_imps_48.pdf');

%performance({'C'}, 1 , [6 12 24 48 96],1e5,1000);
%exportgraphics(gca,'../outputs/perftest_C_6_96.pdf');

%performance(@factory, 1,[6 12 24 48 96],1e5,100);
%exportgraphics(gca,'../outputs/stress_6_96.pdf');

%return

clock_offsets(@factory, 6)
if exp; exportgraphics(gca,'..\outputs\clock_offsets.pdf'); end;

rotation(@factory, 5,2)
if exp; exportgraphics(gca,'..\outputs\rotation2.pdf'); end;
rotation(@factory, 5,1)
if exp; exportgraphics(gca,'..\outputs\rotation1.pdf'); end;

constant(@factory, 1,101,true,false,[NaN,1])
if exp; exportgraphics(gca,'..\outputs\constant_smoothed.pdf'); end;
constant(@factory, 1,101,false,false,[NaN,1])
if exp; exportgraphics(gca,'..\outputs\constant_filtered.pdf'); end;

constant(@factory, 1,101,true,true,[50,0.25])
if exp; exportgraphics(gca,'..\outputs\sloping_smoothed_exception.pdf'); end;
constant(@factory, 1,101,false,true,[50,0.25])
if exp; exportgraphics(gca,'..\outputs\sloping_filtered_exception.pdf'); end;

% projectile
% this function saves the plots to pdf files, so no
% need to call exportgraphics here.
projectile(@factory,1);

% The next script should produce two steps with an estimate near 1,
% than 2 more with one estimate near 1 and the other near 2, and then
% the first parameter is removed, so values near 1 should be returned
% twice more.
add_remove(@factory,1);

function kalman = factory()
    switch implementation
        case 'N'
            kalman = UltimateKalmanNative;
        case 'J'
            kalman = UltimateKalmanJava;
        case 'M'
            kalman = UltimateKalman;
        otherwise
            error("implementation should be 'M', 'N', or 'J' (for MATLAB, native, or Java)");
    end
end

end