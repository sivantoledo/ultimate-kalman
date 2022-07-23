function replication(native)
% test-replication script for UltimateKalman
%
% copyright 2022 Sivan Toledo

if nargin<1
    native = false;
end

addpath 'C:\Users\stoledo\git\ultimate-kalman\native'

%performance(@factory, 1,[6 12],1e6,10000);
%exportgraphics(gca,'../outputs/stress_6_12_long.pdf');

%performance(@factory, 1,[6 12 24 48 96],1e5,1000);
%exportgraphics(gca,'../outputs/stress_6_96.pdf');

%return

clock_offsets(@factory, 6)
exportgraphics(gca,'..\outputs\clock_offsets.pdf');

rotation(@factory, 5,2)
exportgraphics(gca,'..\outputs\rotation2.pdf');
rotation(@factory, 5,1)
exportgraphics(gca,'..\outputs\rotation1.pdf');

constant(@factory, 1,101,true,false,[NaN,1])
exportgraphics(gca,'..\outputs\constant_smoothed.pdf');
constant(@factory, 1,101,false,false,[NaN,1])
exportgraphics(gca,'..\outputs\constant_filtered.pdf');

constant(@factory, 1,101,true,true,[50,0.25])
exportgraphics(gca,'..\outputs\sloping_smoothed_exception.pdf');
constant(@factory, 1,101,false,true,[50,0.25])
exportgraphics(gca,'..\outputs\sloping_filtered_exception.pdf');

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
    if native
        kalman = UltimateKalmanNative;
    else
        kalman = UltimateKalman;
    end
end

end