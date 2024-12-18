function [initConds,dataMats,sysFuncs] = initialStructs(j,k,m,n,QRH)
%initialStructs Takes in the number of observed and state variables, the number of
% satellites and timesteps as well as whether or not to preallocate various 
% filter matrices and outputs appropriately sized structures with the necessary
% fields to be set and facilitate simulation.
% INPUTS:
%      j - Number of observed variables
%      k - Number of state variables
%      m - Number of satellites or constellations
%      n - Number of time steps
%    QRH - Structure of Q, R, and H data and types
%        . Qdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Rdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Hdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Qtype - string telling if const, data, SNC, DMC, or other func
%        . Rtype - string telling if const, data, or func
%        . Htype - string telling if const, data, or func
% OUTPUT:
% sConstell - structure of required parameters for generating initial orbits
%           .   satFun - Satellite constellation generation function handle
%           .   radius - Radius of orbits in kilometers
%           .   inclin - Inclination of orbits in degrees
%           .    nSats - Total number of satellites
%           .  nPlanes - Number of geometry planes
%           .  phasing - Phasing between satellites
%           . satNames - String describing name of satellite constellation
% initConds - Structure containing initial conditions
%           .  X_est0 - Initial state estimate
%           . dx_est0 - Initial state uncertainty estimate
%           .     P_0 - Initial covariance estimate
%  dataMats - Structure containing data matrices
%        . Qdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Rdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Hdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Qtype - string telling if const, data, SNC, DMC, or other func
%        . Rtype - string telling if const, data, or func
%        . Htype - string telling if const, data, or func
%  sysFuncs - Structure containing system functions
%           .  XdotPhidot - State derivative function handle
%           .    vsIntFun - Variable step integration function handle
%           .     odeOpts - Default odeset options for vsIntFun
% TODO check that default functions exist
arguments (Input)
    j (1,1) {mustBeNumeric,mustBeNonempty,mustBeFinite}
    k (1,1) {mustBeNumeric,mustBeNonempty,mustBeFinite}
    m (1,1) {mustBeNumeric,mustBeNonempty,mustBeFinite}
    n (1,1) {mustBeNumeric,mustBeNonempty,mustBeFinite}
    QRH (1,1) struct {mustBeValidQRH}
end

% Temp
I = eye(k);
I(I==1) = Inf;

% Initial Conditions
initConds.X_est0 = zeros(m,k);    % Default value is zero, be sure to set
initConds.dx_est0 = zeros(m,k);   % Default value should be zero vector
initConds.P_0 = I + zeros(k,k,m); % Default value should be infinite

% Data Matrices
dataMats.Qtype = QRH.Qtype;
dataMats.Rtype = QRH.Rtype;
dataMats.Htype = QRH.Htype;
if isfield(QRH,"Qdata")
    dataMats.Qdata = QRH.Qdata;
else
    switch QRH.Qtype
        case "const"
            dataMats.Qdata = zeros(k,k,m);
        case "data"
            % Avoid Using Unless your system is relatively simple
            dataMats.Qdata = zeros(k,k,n,m);
        case "SNC"
            dataMats.Qdata = zeros(j,j,m);
        case "DMC"
            dataMats.Qdata = zeros(j,j,m);
        case "func"
            dataMats.Qdata = @Q_matrix_function_handle;
        otherwise
            eid = "QMatrix:invalidQtype";
            msg = "This shouldn't be possible.";
            error(eid,msg)
    end
end
if isfield(QRH,"Rdata")
    dataMats.Rdata = QRH.Rdata;
else
    switch QRH.Rtype
        case "const"
            dataMats.Rdata = zeros(j,j,m);
        case "data"
            % Avoid Using Unless your system is relatively simple
            dataMats.Rdata = zeros(j,j,n,m);
        case "func"
            dataMats.Rdata = @R_matrix_function_handle;
        otherwise
            eid = "RMatrix:invalidRtype";
            msg = "This shouldn't be possible.";
            error(eid,msg)
    end
end
if isfield(QRH,"Hdata")
    dataMats.Hdata = QRH.Hdata;
else
    switch QRH.Htype
        case "const"
            dataMats.Hdata = zeros(j,k,m);
        case "data"
            % Avoid Using Unless your system is relatively simple
            dataMats.Hdata = zeros(j,k,n,m);
        case "func"
            dataMats.Hdata = @H_matrix_function_handle;
        otherwise
            eid = "HMatrix:invalidHtype";
            msg = "This shouldn't be possible.";
            error(eid,msg)
    end
end

% System Functions
sysFuncs.XdotPhidot = @stateSTM_nodrag;                 % Default value is simple, no drag
sysFuncs.vsIntFun = @ode45;                             % Default value is ode45
sysFuncs.odeOpts = odeset("AbsTol",1e-9,"RelTol",1e-8); % Default tolerances

end % initialStructs Function

%% Validation Functions
function mustBe2Dor3DArrayorFunc(A)
    % Checks if A is a 2D or 3D array or a function
    if isa(A,"function_handle")
        if ~(isfile(func2str(A) + ".m") || ...
            exist(func2str(A),"file") == 2 || ...
            exist(func2str(A),"builtin") == 5)
            eid = "Type:notAValidFunctionHandle";
            msg = "Must be a valid function handle. Check for typos and " + ...
            "ensure that the file/function is on your path.";
            error(eid,msg)
        end
    elseif isa(A,"double")
        if ~((ismatrix(A) || ndims(A) == 3) && isa(A,"double"))
            eid = "Type:not2Dor3DArray";
            msg = "Must be a 2D or 3D array of doubles.";
            error(eid,msg)
        end
    else
        eid = "Type:not2Dor3DArrayorFunction";
        msg = "Must be a 2D or 3D array of doubles or a function handle.";
        error(eid,msg)
    end
end

function mustBeValidQRH(QRH)
    % Checks to see if the fields of QRH are correct
    if isfield(QRH,"Qdata")
        mustBe2Dor3DArrayorFunc(QRH.Qdata)
    end
    if isfield(QRH,"Rdata")
        mustBe2Dor3DArrayorFunc(QRH.Rdata)
    end
    if isfield(QRH,"Hdata")
        mustBe2Dor3DArrayorFunc(QRH.Hdata)
    end
    mustBeTextScalar(QRH.Qtype)
    mustBeMember(QRH.Qtype,["const","data","SNC","DMC","func"])
    mustBeTextScalar(QRH.Rtype)
    mustBeMember(QRH.Rtype,["const","data","func"])
    mustBeTextScalar(QRH.Htype)
    mustBeMember(QRH.Htype,["const","data","func"])
end