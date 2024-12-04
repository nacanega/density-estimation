function sConstell = initialConstell(constName)
%initialConstell returns a default satellite constellation structure
% INPUT:
% constname - String (optional) Defaults to just initialized structure
%             If specified, will load an existing constellation if a match
%             Iridium, Galileo, and generic Walker Star/Delta
% OUTPUT:
% sConstell - structure of required parameters for generating initial orbits
%           .   satFun - Satellite constellation generation function handle
%           .   radius - Radius of orbits in kilometers
%           .   inclin - Inclination of orbits in degrees
%           .    nSats - Total number of satellites
%           .  nPlanes - Number of geometry planes
%           .  phasing - Phasing between satellites
%           .   argLat - Argument of latitude in degrees
%           . satNames - String describing name of satellite constellation
%         - Alternate arguments for single satellite
%           .   satFun - Satellite constellation generation function handle
%           .semimajor - Semimajor axis or orbit [km]
%           .   eccent - Eccentricity
%           .   inclin - Inclination of orbits in [deg]
%           .     RAAN - Right Ascension of ascending node [deg]
%           .  argPeri - Argument of periapsis [deg]
%           . trueAnom - True Anomaly [deg]
%           .    nSats - Total number of satellites
%           . satNames - String describing name of satellite constellation
arguments (Input)
    constName string = "Default"
end

switch constName
    case {"single","Single"}
        % Single Satellite
        sConstell.satFun = @satellite;    % Single satellite orbit
        sConstell.semimajor = 0;          % Radius in [km]
        sConstell.eccent = 0;             % Eccentricity
        sConstell.inclin = 0;             % Inclination [deg]
        sConstell.RAAN = 0;               % Right ascension of ascending node [deg]
        sConstell.argPeri = 0;            % Argument of periapsis [deg]
        sConstell.trueAnom = 0;           % True anomaly [deg]
        sConstell.nSats = 1;              % Number of satellites
        sConstell.satNames = "Satellite"; % Satellite name
    case {"singleCircular","SingleCircular"}
        % Single Satellite Circular
        sConstell.satFun = @satellite;   % Single satellite orbit
        sConstell.semimajor = 7159.1363; % Radius in [km]
        sConstell.eccent = 0;            % Eccentricity
        sConstell.inclin = 28.5;         % Inclination [deg]
        sConstell.RAAN = 0;              % Right ascension of ascending node [deg]
        sConstell.argPeri = 0;           % Argument of periapsis [deg]
        sConstell.trueAnom = 0;          % True anomaly [deg]
        sConstell.nSats = 1;             % Number of satellites
        sConstell.satNames = "Circular"; % Satellite name
    case {"singleElliptical","SingleElliptical"}
        % Single Satellite Elliptical
        sConstell.satFun = @satellite;     % Single satellite orbit
        sConstell.semimajor = 8378.1363;   % Radius in [km]
        sConstell.eccent = 0.2;            % Eccentricity
        sConstell.inclin = 28.5;           % Inclination [deg]
        sConstell.RAAN = 0;                % Right ascension of ascending node [deg]
        sConstell.argPeri = 0;             % Argument of periapsis [deg]
        sConstell.trueAnom = 0;            % True anomaly [deg]
        sConstell.nSats = 1;               % Number of satellites
        sConstell.satNames = "Elliptical"; % Satellite name
    case {"doubleCircular","DoubleCircular"}
        % Double Satellite Circular
        sConstell.satFun = @nSat;        % N satellite orbit
        sConstell.semimajor = 7159.1363; % Radius in [km]
        sConstell.eccent = 0;            % Eccentricity
        sConstell.inclin = 28.5;         % Inclination [deg]
        sConstell.RAAN = [0;180];        % Right ascension of ascending node [deg]
        sConstell.argPeri = 0;           % Argument of periapsis [deg]
        sConstell.trueAnom = 0;          % True anomaly [deg]
        sConstell.nSats = 2;             % Number of satellites
        sConstell.satNames = "Circular"; % Satellite name
    case {"nCircular","NCircular"}
        % N Satellite Circular
        sConstell.satFun = @nSat;        % N satellite orbit
        sConstell.semimajor = 7159.1363; % Radius in [km]
        sConstell.eccent = 0;            % Eccentricity
        sConstell.inclin = 28.5;         % Inclination [deg]
        sConstell.RAAN = 0;              % Right ascension of ascending node [deg]
        sConstell.argPeri = 0;           % Argument of periapsis [deg]
        sConstell.trueAnom = 0;          % True anomaly [deg]
        sConstell.nSats = 0;             % Number of satellites
        sConstell.satNames = "Circular"; % Satellite name
    case {"Iridium","iridium"}
        % Satellite Constellation (Iridium)
        sConstell.satFun = @walkerStar; % Iridium has a Walker Star configuration
        sConstell.radius = 7159.1363;   % Radius in [km]
        sConstell.inclin = 86.4;        % Inclination [deg]
        sConstell.nSats = 66;           % Number of satellites
        sConstell.nPlanes = 6;          % Number of planes
        sConstell.phasing = 2;          % Phasing
        sConstell.argLat = 0;           % Argument of latitude [deg]
        sConstell.satNames = "Iridium"; % Name String
    case {"Galileo","galileo"}
        % Satellite Constellation (Galileo)
        sConstell.satFun = @walkerDelta; % Galileo has a Walker Delta configuration
        sConstell.radius = 29599.8;      % Radius in [km]
        sConstell.inclin = 56;           % Inclination [deg]
        sConstell.nSats = 24;            % Number of satellites
        sConstell.nPlanes = 3;           % Number of planes
        sConstell.phasing = 1;           % Phasing
        sConstell.argLat = 15;           % Argument of latitude [deg]
        sConstell.satNames = "Galileo";  % Name String
    case {"walkerStar","star","Star","ws"}
        % Satellite Constellation (Walker Star)
        sConstell.satFun = @walkerStar; % Walker Star configuration
        sConstell.radius = 0;           % Radius in [km]
        sConstell.inclin = 0;           % Inclination [deg]
        sConstell.nSats = 0;            % Number of satellites
        sConstell.nPlanes = 0;          % Number of planes
        sConstell.phasing = 0;          % Phasing
        sConstell.argLat = 0;           % Argument of latitude [deg]
        sConstell.satNames = "satName"; % Name String
    case {"walkerDelta","delta","Delta","wd"}
        % Satellite Constellation (Walker Delta)
        sConstell.satFun = @walkerDelta; % Walker Delta configuration
        sConstell.radius = 0;           % Radius in [km]
        sConstell.inclin = 0;           % Inclination [deg]
        sConstell.nSats = 0;            % Number of satellites
        sConstell.nPlanes = 0;          % Number of planes
        sConstell.phasing = 0;          % Phasing
        sConstell.argLat = 0;           % Argument of latitude [deg]
        sConstell.satNames = "satName"; % Name String
    otherwise
        % Satellite Constellation (Default)
        sConstell.satFun = @customFunc; % Custom configuration function
        sConstell.radius = 0;           % Radius in [km]
        sConstell.inclin = 0;           % Inclination [deg]
        sConstell.nSats = 0;            % Number of satellites
        sConstell.nPlanes = 0;          % Number of planes
        sConstell.phasing = 0;          % Phasing
        sConstell.argLat = 0;           % Argument of latitude [deg]
        sConstell.satNames = "satName"; % Name String
end

end