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
arguments (Input)
    constName string = "Default"
end

switch constName
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
        % Satellite Constellation (Iridium)
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