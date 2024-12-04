function [satScen,initElems,initStates,varargout] = initialOrbits(S,satStates,modStates,satParams,defParams)
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
% TODO add validation
arguments (Input)
    S struct
    satStates string {mustBeVector} = ["r_x","r_y","r_z","v_x","v_y","v_z"];
    modStates string {mustBeVector} = "";
    satParams string {mustBeVector} = "";
    defParams double = [];
end
narginchk(1,5)
nargoutchk(3,5)

% startDate = datetime(2000,1,1,12,0,0);
% stopDate = startDate + hours(1.5);
% sampleTime = 60;
satScen = satelliteScenario("AutoSimulate",false);
satFun = S.satFun;
satFunName = func2str(satFun);

switch satFunName
    case {"walkerDelta", "walkerStar"}
        % walkerDelta or walkerStar
        sats = satFun(satScen, 1e3*S.radius, S.inclin, S.nSats, S.nPlanes, ...
            S.phasing, ArgumentofLatitude=S.argLat, Name=S.satNames);
    case "satellite"
        % singleSat
        sats = satFun(satScen, 1e3*S.semimajor, S.eccent, S.inclin,...
            S.RAAN, S.argPeri, S.trueAnom, Name=S.satNames);
    otherwise
        % Custom satFun
        sats = satFun(satScen,S);
end

% Iterate over each satellite and obtain elements
nSats = S.nSats;
advance(satScen);

if ~(satParams == "")
    nParams = length(satParams);
    initParams = ones(nSats,nParams);
    if ~isempty(defParams)
        [nSatChk,~] = size(defParams);
        if nSatChk == 1
            defParams = defParams.*ones(nSats,nParams);
        elseif nSatChk ~= nSats
            eid = "Size:specifiedParametersMustBeVectorOrMatchSatelliteNumber";
            msg = "Defined parameters must either be a single vector or " + ...
                "a matrix where the number of rows matches the satellite number.";
            error(eid,msg);
        end
    end
else
    nParams = 0;
end

for i = nSats:-1:1
    % Set orbital elements
    OEs = orbitalElements(sats(i));

    % Common Elements
    initElems(i,7) = OEs.Period;
    initElems(i,2) = OEs.Eccentricity;
    initElems(i,3) = OEs.Inclination;
    initElems(i,4) = OEs.RightAscensionOfAscendingNode;
    initElems(i,5) = OEs.ArgumentOfPeriapsis;
    
    % Different based on number of satellites...
    if isfield(OEs,"TrueAnomaly")
        initElems(i,1) = OEs.SemiMajorAxis;
        initElems(i,6) = OEs.TrueAnomaly;
    else
        initElems(i,1) = OEs.MeanMotion;
        initElems(i,6) = OEs.MeanAnomaly;
    end

    % Set states
    switch satStates
        case ["r_x","r_y","r_z","v_x","v_y","v_z"]
            [pos, vel] = states(sats(i));
            pos = pos*1e-3; vel = vel*1e-3;
            initStates(i,:) = [pos(:,1).' vel(:,1).'];
        case ["r_x","r_y","r_z","v_x","v_y","v_z","C_D"]
            [pos, vel] = states(sats(i));
            pos = pos*1e-3; vel = vel*1e-3;
            initStates(i,:) = [pos(:,1).' vel(:,1).' ...
                sats(i).PhysicalProperties.DragCoefficient];
        otherwise
            eid = "States:undefinedStateSequence";
            msg = "State vector not recognized, please check for typos, " + ...
                "edit this file to include the case. Alternatively, " + ...
                "construct your initial states manually by calling " + ...
                "without the state, model, and parameter vectors.";
            error(eid,msg)
    end

    if nargout >= 4
        % Set model states
        % TODO Nonspherical Earth
        % TODO Circular and elliptical orbit state native
        h = norm(pos) - 6378.1363;
        switch modStates
            case ""
                modVals = [];
            case "rho_0" % Exponential model single parameter estimation
                [rho_0,h_0,H] = expParams(h);
                modVals(i) = rho_0;
            case ["rho_0","H"] % Requires satellites at different altitudes
                [rho_0,h_0,H] = expParams(h);
                modVals(i,2) = H;
                modVals(i,1) = rho_0;
            case ["rho_0","h_0","H"] % Requires satellites at different altitudes
                [rho_0,h_0,H] = expParams(h);
                modVals(i,3) = H;
                modVals(i,2) = h_0;
                modVals(i,1) = rho_0;
            otherwise
                eid = "Model:undefinedStateSequence";
                msg = "Model vector not recognized, please check for typos, " + ...
                    "edit this file to include the case. Alternatively, " + ...
                    "construct your initial states manually by calling " + ...
                    "without the state, model, and parameter vectors.";
                error(eid,msg)
        end
        if nargout == 5
            % Parameter output
            if nParams > 0 && isempty(defParams)
                % Parameters are specified, but not defined
                for j = length(satParams):-1:1
                    pName = satParams(j);
                    switch pName
                        case {"A","A_drag","A_DRAG"}
                            initParams(i,j) = ... 
                                sats(i).PhysicalProperties.DragArea;
                        case {"C_D","C_d"}
                            initParams(i,j) = ...
                                sats(i).PhysicalProperties.DragCoefficient;
                        case {"A_SRP","A_srp"}
                            initParams(i,j) = ...
                                sats(i).PhysicalProperties.SRPArea;
                        case {"C_R","C_r"}
                            initParams(i,j) = ...
                                sats(i).PhysicalProperties.ReflectivityCoefficient;
                        case {"m","m_s","mass","Mass","MASS"}
                            initParams(i,j) = ...
                                sats(i).PhysicalProperties.Mass;
                        case "h_0"
                            initParams(i,j) = h_0;
                        case "H"
                            initParams(i,j) = H;
                        otherwise
                            eid = "Parameter:undefinedSatelliteParameter";
                            error(eid, ...
                                "'%s' is an undefined satellite parameter.",pName)
                    end
                end
            elseif nParams > 0 && ~isempty(defParams)
                % Use specified parameters
                initParams(i,:) = defParams(i,:);
            else % No parameters are specified or defined, use matlab defaults
                initParams(i,5) = sats(i).PhysicalProperties.Mass;
                initParams(i,4) = sats(i).PhysicalProperties.SRPArea;
                initParams(i,3) = sats(i).PhysicalProperties.ReflectivityCoefficient;
                initParams(i,2) = sats(i).PhysicalProperties.DragArea;
                initParams(i,1) = sats(i).PhysicalProperties.DragCoefficient;
            end
        end
    end
end % for

if nargout == 5
    varargout{2} = initParams;
    varargout{1} = modVals;
elseif nargout == 4
    varargout{1} = modVals;
end

end % function