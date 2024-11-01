function [satScen,initElems,initStates] = initialOrbits(S,varargin)
%initialOrbits takes in an initial orbit generation structure consisting of
% the function name, the required arguments, and names for the satellites and
% outputs a handle to the satellite scenario, initial orbital elements, the
% initial state
% INPUTS:
% S - structure containing required parameters for generating initial orbits
%     .   satFun - Satellite constellation generation function handle
%     .   radius - Radius of orbits in kilometers
%     .   inclin - Inclination of orbits in degrees
%     .    nSats - Total number of satellites
%     .  nPlanes - Number of geometry planes
%     .  phasing - Phasing between satellites
%     .   argLat - Argument of latitude in degrees
%     . satNames - String describing name of satellite constellation
% varargin - {1} [1,m] String vector of state variable names
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
% initStates - [n,m] Matrix of initial states for each satellite
%
% note: this will change when drag is added, have a separate function for drag
% TODO Account for drag and density parameters

satScen = satelliteScenario;
satFun = S.satFun;
satFunName = func2str(satFun);

if strcmp(satFunName,"walkerDelta") || strcmp(satFunName,"walkerStar")
    % walkerDelta or walkerStar
    sats = satFun(satScen, 1e3*S.radius, S.inclin, S.nSats, S.nPlanes, ...
        S.phasing, ArgumentofLatitude=S.argLat, Name=S.satNames);
else
    % Custom satFun
    sats = satFun(satScen,S);
end

% Iterate over each satellite and obtain elements
nSats = length(sats);

for i = nSats:-1:1
    OEs = orbitalElements(sats(i));
    initElems(i,1) = OEs.SemiMajorAxis;
    initElems(i,2) = OEs.Eccentricity;
    initElems(i,3) = OEs.Inclination;
    initElems(i,4) = OEs.RightAscensionOfAscendingNode;
    initElems(i,5) = OEs.ArgumentOfPeriapsis;
    initElems(i,6) = OEs.TrueAnomaly;
    initElems(i,7) = OEs.Period;
    [pos, vel] = states(sats(i));
    initStates(i,:) = [pos(:,1).' vel(:,1).']./1e3;
end

end