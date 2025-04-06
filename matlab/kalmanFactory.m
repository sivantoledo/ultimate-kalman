function factory = kalmanFactory(className, options)
%KALMANFACTORY Create a function handle to instantiate a Kalman filter class object.
%
%   FACTORY = KALMANFACTORY(CLASSNAME, OPTIONS) returns a function handle that,
%   when called with no arguments, instantiates an object of the specified Kalman
%   filter handle class using the provided options struct. This function serves
%   as a factory for creating instances of Kalman filter classes, such as
%   'KalmanAssociativeSmoother', 'KalmanConventional', etc., each of which must
%   be a handle class with a constructor that accepts a single struct argument.
%
%   FACTORY = KALMANFACTORY(CLASSNAME) uses an empty struct as the default options.
%
%   Inputs:
%       CLASSNAME - String scalar specifying the name of the Kalman filter class
%                   to instantiate. Must be one of the following (case-sensitive):
%                   'KalmanAssociativeSmoother', 'KalmanConventional',
%                   'KalmanExplicitRepresentation', 'KalmanJava', 'KalmanNative',
%                   'KalmanOddevenSmoother', 'KalmanSparse', 'KalmanUltimate'.
%       OPTIONS   - Struct containing options to pass to the class constructor.
%                   Optional; defaults to an empty struct if not provided. The
%                   fields of OPTIONS must match the expectations of the specified
%                   class's constructor.
%
%   Output:
%       FACTORY - Function handle that, when called as FACTORY(), returns a new
%                 instance of the specified Kalman filter class initialized with
%                 OPTIONS.
if nargin<2
    options = struct();
end

kalmanClasses = {'KalmanAssociativeSmoother', 'KalmanConventional', 'KalmanPython', 'KalmanExplicitRepresentation', 'KalmanJava', 'KalmanNative', 'KalmanOddevenSmoother', 'KalmanSparse', 'KalmanUltimate'}';

if ~ischar(className) || ~ismember(className, kalmanClasses)
    names='';
    for i=1:length(kalmanClasses)
        names=[names ' ' char(kalmanClasses{i})];
    end
    error('className must be one of: %s',names);
end
if ~isstruct(options)
    error('params must be a struct');
end

factory = @() feval(className, options);

end

