function [satScen,initElems,initStates,varargout] = initialStates(S,stateType,satStates,modStates,satParams,defParams)
%initialOrbits takes in an initial orbit generation structure consisting of
% the function name, the required arguments, and names for the satellites and
% outputs a handle to the satellite scenario, initial orbital elements, the
% initial state
% INPUTS:
% S - structure containing required parameters for generating initial orbits
%     .   satFun - Satellite constellation generation function handle
%     .   radius - Radius of orbits in [km]
%     .   inclin - Inclination of orbits in [deg]
%     .    nSats - Total number of satellites
%     .  nPlanes - Number of geometry planes
%     .  phasing - Phasing between satellites
%     .   argLat - Argument of latitude in degrees
%     . satNames - String describing name of satellite constellation
%   - Alternate arguments for single satellite
%     .   satFun - Satellite constellation generation function handle
%     .semimajor - Semimajor axis or orbit [km]
%     .   eccent - Eccentricity
%     .   inclin - Inclination of orbits in [deg]
%     .     RAAN - Right Ascension of ascending node [deg]
%     .  argPeri - Argument of periapsis [deg]
%     . trueAnom - True Anomaly [deg]
%     .    nSats - Total number of satellites
%     . satNames - String describing name of satellite constellation
% stateType - String of type of states to output
% satStates - (1,m) String vector of state variable names
% modStates - (1,k) String vector of model state variable names
% satParams - (1,p) String vector of system parameters
% defParams - (1,p) Vector of default values for system parameters
% OUTPUT:
%    satScen - Satellite scenario object
%  initElems - [n,7] Matrix of initial orbital elements with the columns:
%              [a, e, i, RAAN, omega, nu, P]
%                 a - Semimajor axis
%                 e - Eccentricity
%                 i - Inclination
%              RAAN - Right ascension of the ascending node
%             omega - Argument of periapsis
%                nu - True anomaly
%                 P - Period
% SingleSat - [n, e, i, RAAN, omega, M, P]
%                 n - Mean motion
%                 e - Eccentricity
%                 i - Inclination
%              RAAN - Right ascension of the ascending node
%             omega - Argument of periapsis
%                 M - Mean anomaly
%                 P - Period
% initStates - [n,m] Matrix of initial states for each satellite
%  varargout - {1} [n,p] Matrix of parameters for each satellite
% Note that (n) may equal 1 if all satellites are included in the same
% state vector such as when estimating parameters
% TODO add validation
arguments (Input)
    S struct
    stateType (1,1) string = "separate"; 
    satStates string {mustBeVector} = ["r_x","r_y","r_z","v_x","v_y","v_z"];
    modStates string {mustBeVector} = "";
    satParams string {mustBeVector} = "";
    defParams double = [];
end

% Setup input arguments
if ~(defParams == "")
    if ~(satParams == "")
        initArgs = {satStates, modStates, satParams, defParams};
    end
else
    if ~(satParams == "")
        initArgs = {satStates, modStates, satParams};
    end
end

if stateType == "aggregate" || stateType == "combined"
    % We will need to combine the obtained data to a single vectors

    % Obtain separate satellite data
    if nargout == 4
        [satScen,sepElems,sepStates,sepMod,varargout{1}] = initialOrbits(S,initArgs{:});
    else
        [satScen,sepElems,sepStates,sepMod] = initialOrbits(S,initArgs{:});
    end % if

    m = S.nSats;
    n = length(satStates);
    o = length(modStates);

    initElems = zeros(6*m + 1,1);
    initStates = zeros(n*m + o,1);

    for i = 1:m
        initElems(6*(i-1)+1:6*i) = sepElems(i,1:6)';
        initStates(n*(i-1)+1:n*i) = sepStates(i,:)';
    end

    initElems(end) = mean(sepElems(:,7),"all"); % Average Orbital Period
    initStates(m+1:end) = mean(sepMod,1)'; % Average Model parameters, assume all sats share
    
else % Individual satellites

    % Obtain separate satellite data
    if nargout == 4
        [satScen,initElems,sepStates,sepMod,varargout{1}] = initialOrbits(S,initArgs{:});
    else
        [satScen,initElems,sepStates,sepMod] = initialOrbits(S,initArgs{:});
    end % if

    initStates = [sepStates, sepMod];

end % if

end % function