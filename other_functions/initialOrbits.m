function [satScen,initElems,initStates,varargout] = initialOrbits(S,varargin)
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
%          - {2} [1,p] String vector of system parameters
%          - {3} [1,p] Vector of default values for system parameters
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
%  varargout - {1} [n,p] Matrix of parameters for each satellite
% TODO add validation
narginchk(1,4)
nargoutchk(3,4)

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

% Set orbital elements
for i = nSats:-1:1
    OEs = orbitalElements(sats(i));
    initElems(i,1) = OEs.SemiMajorAxis;
    initElems(i,2) = OEs.Eccentricity;
    initElems(i,3) = OEs.Inclination;
    initElems(i,4) = OEs.RightAscensionOfAscendingNode;
    initElems(i,5) = OEs.ArgumentOfPeriapsis;
    initElems(i,6) = OEs.TrueAnomaly;
    initElems(i,7) = OEs.Period;
end

% Different states depending on how many states we have
if nargin > 1
    % We have a different number of states defined in varagin{1}
    switch varargin{1}
        case ["r_x","r_y","r_z","v_x","v_y","v_z"]
            % We can use the standard states
            for i = nSats:-1:1
                [pos, vel] = states(sats(i));
                initStates(i,:) = [pos(:,1).' vel(:,1).']./1e3;
            end
        case ["r_x","r_y","r_z","v_x","v_y","v_z","rho_0","h_0","H"]
            for i = nSats:-1:1
                [pos, vel] = states(sats(i));
                pos = pos*1e-3; vel = vel*1e-3;
                r = norm(pos(:,1));
                modPs = densityParams(r);
                initStates(i,:) = [pos(:,1).' vel(:,1).' modPs];
            end
        otherwise
            eid = "States:undefinedStateSequence";
            msg = "State vector not recognized, please check for typos, " + ...
                "edit this file to include the case. Alternatively, " + ...
                "construct your initial states manually by calling " + ...
                "without the state and parameter vectors.";
            error(eid,msg)
    end
else
    % We can use the standard states
    for i = nSats:-1:1
        [pos, vel] = states(sats(i));
        initStates(i,:) = [pos(:,1).' vel(:,1).']./1e3;
    end
end

% Different parameter outputs depending on which are specified
if nargout == 4
    % We are outputting parameters
    if nargin == 1 || nargin == 2
        % Parameters are not specified, output default
        % Cd A m 
        for i = nSats:-1:1
            initParams(i,5) = sats(i).PhysicalProperties.Mass;
            initParams(i,4) = sats(i).PhysicalProperties.SRPArea;
            initParams(i,3) = sats(i).PhysicalProperties.ReflectivityCoefficient;
            initParams(i,2) = sats(i).PhysicalProperties.DragArea;
            initParams(i,1) = sats(i).PhysicalProperties.DragCoefficient;
        end
    else
        % Parameters are specified
        nParams = length(varargin{2});
        if nargin == 4
            % Use specified parameters
            initParams = ones(nSats,nParams);
            initParams = initParams.*varargin{3};
        else
            % Use default paramters
            for i = nSats:-1:1
                for j = nParams:-1:1
                    switch varargin{2}(j)
                        case {"C_D","CD","DragCoefficient"}
                            pName = "DragCoefficient";
                        case {"A","A_D","AD","Area","area","DragArea"}
                            pName = "DragArea";
                        case {"C_R","CR","ReflectivityCoefficient"}
                            pName = "ReflectivityCoefficient";
                        case {"A_R","AR","SRPArea"}
                            pName = "SRPArea";
                        case {"m","m_s","m_S","M","M_s","M_S","mass","Mass"}
                            pName = "Mass";
                        otherwise
                            eid = "Parameter:undefinedSatelliteParameter";
                            error(eid, ...
                                "'%s' is an undefined satellite parameter", ...
                                varargin{2}(j))
                    end
                    initParams(i,j) = sats(i).PhysicalProperties.(pName);
                end
            end
        end
    end
    varargout{1} = initParams;
end