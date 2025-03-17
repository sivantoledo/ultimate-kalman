function kalman = kalmanFactory(class,options)

if nargin<2
    options = struct();
end

switch class
    case 'KalmanNative'
        kalman = KalmanNative(options);
    case 'KalmanJava'
        kalman = KalmanJava(options);
    case 'KalmanUltimate'
        kalman = KalmanUltimate(options);
    case 'KalmanConventional'
        kalman = KalmanConventional(options);
    case 'KalmanAssociativeSmoother'
        kalman = KalmanAssociativeSmoother(options);
    case 'KalmanOddevenSmoother'
        kalman = KalmanOddevenSmoother(options);
    case 'KalmanSparse'
        kalman = KalmanSparse(options);
    otherwise
        error(sprintf('%s is not a Kalman class',class));
end


end

