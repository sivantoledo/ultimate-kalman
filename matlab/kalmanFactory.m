function kalman = kalmanFactory(class,options)

if nargin<2
    options = struct();
end

switch class
    case 'UltimateKalmanNative'
        kalman = UltimateKalmanNative(options);
    case 'UltimateKalmanJava'
        kalman = UltimateKalmanJava(options);
    case 'UltimateKalman'
        kalman = UltimateKalman(options);
    case 'KalmanFilterSmoother'
        kalman = KalmanFilterSmoother(options);
    case 'KalmanAssociativeSmoother'
        kalman = KalmanAssociativeSmoother(options);
    case 'KalmanSparse'
        kalman = KalmanSparse(options);
    otherwise
        error(sprintf('%s is not a Kalman class',class));
end


end

