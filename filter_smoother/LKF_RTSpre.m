function [filSol,smoSol,addOut] = LKF_RTSpre(t,zs,initConds,dataMats,sysFuncs,opts)
%LKF_RTSpre Is a more costly version of an LKF_RTS filter and smoother
% which precalculates most of the matrices outside of the main filter.
% INPUTS:
%         t - [n,1] Vector of times
%        zs - [n,k] Matrix of observations
% initConds - Structure if initial conditions
%           .  X_est0 - Initial state estimate
%           . dx_est0 - Initial state uncertainty estimate
%           .     P_0 - Initial covariance estimate
% TODO      .  params - Input for state functions (ie drag, area, etc)
%  dataMats - Structure of data matrices
%        . Qdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Rdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Hdata - 2D or 3D matrix of doubles or function handle (optional)
%        . Qtype - string telling if const, data, SNC, DMC, or other func
%        . Rtype - string telling if const, data, or func
%        . Htype - string telling if const, data, or func
%  sysFuncs - Structure of system functions
%           .  XdotPhidot - State derivative function handle
%           .    vsIntFun - Variable step integration function handle
%           .     odeOpts - Default odeset options for vsIntFun
%           .      params - Structure of parameters for state equations
%      opts - Structure of filter options
%           .     tol - tolerance to stop iteratively smoothing
%           . maxIter - maximum number of iterations for filter/smoother
%           .  maxInc - maximum number of consecutive increasing iterations
%           . outIter - "last" or "all", which iterations to output
%           . outPmat - "none" or "all", which covariance matrices to output
% OUTPUT:
% filSol - Structure containing the linearized Kalman filter solution
%        .  X - The filtered/estimated state
%        . dx - The filtered/estimated state uncertainty
%        . bs - The filter innovation vectors
%        .  P - The filtered/estimated covariance
% smoSol - Structure containting the iteratively smoothed solution
%        .  X - The smoothed state
%        . dx - The smoothed uncertainty
%        .  P - The smoothed covariance
% addOut - Structure of additional outputs
%        .      iter - The number of iterations performed
%        .   maxDiff - The difference metrix for the initial state estimate
%        . exitState - Why the filter stopped iterating

% FIXME Add parameters input validation

% Input Argument Validation
arguments (Input)
    t double {mustBeNonempty,mustBeFinite,mustBeColumn}
    zs double {mustBeMatrix,observationsMustMatchTimeLength(t,zs)}
    initConds (1,1) struct {mustBeValidInitConds(initConds)}
    dataMats (1,1) struct {mustBeValidDataMats(dataMats)}
    sysFuncs (1,1) struct {mustBeValidSysFuncs(sysFuncs)}
    opts (1,1) struct {mustBeValidOpts(opts)}
end

% Load initial conditions
X_est0 = initConds.X_est0;
dx_est0 = initConds.dx_est0;
P_0 = initConds.P_0;
% Swtiched to normalized maximum difference
order = 10.^floor(log10(X_est0));

% Get sizes for indexing and preallocation
M = size(X_est0,1);
[N,K] = size(zs);

% Load data matrices (or preallocate)
Qtype = dataMats.Qtype;
Rtype = dataMats.Rtype;
Htype = dataMats.Htype;
Qdata = dataMats.Qdata;
Rdata = dataMats.Rdata;
Hdata = dataMats.Hdata;

switch Qtype
    case {"const","data"}
        % Load Q data
        Qs = Qdata;
        Qfunc = false;
    case {"SNC","DMC","func"}
        % Preallocate Q
        Qs = zeros(M,M,N);
        Qfunc = true;
    otherwise
        eid = "QMatrix:invalidQtype";
        msg = "This shouldn't be possible.";
        error(eid,msg)
end

switch Rtype
    case {"const","data"}
        % Load R data
        Rs = Rdata;
        Rfunc = false;
    case "func"
        % Preallocate R
        Rs = zeros(K,K,N);
        Rfunc = true;
    otherwise
        eid = "RMatrix:invalidRtype";
        msg = "This shouldn't be possible.";
        error(eid,msg)
end

switch Htype
    case {"const","data"}
        % Load H data
        Hs = Hdata;
        Hfunc = false;
    case "func"
        % Preallocate H
        Hs = zeros(K,M,N);
        Hfunc = true;
    otherwise
        eid = "HMatrix:invalidHtype";
        msg = "This shouldn't be possible.";
        error(eid,msg)
end   
    
% Load system functions
XdotPhidot = sysFuncs.XdotPhidot;
vsIntFun = sysFuncs.vsIntFun;
odeOpts = sysFuncs.odeOpts;

if nargin(XdotPhidot) == 2
    dotFunc = XdotPhidot;
elseif nargin(XdotPhidot) == 3
    dotFunc = @(t,X) XdotPhidot(t,X,sysFuncs.params); 
else
    eid = "Input:invalidArgumentNumber";
    msg = "This shouldn't be possible.";
    error(eid,msg);
end

% Load Filter Options
tol = opts.tol;
maxIter = opts.maxIter;
maxInc = opts.maxInc;
outIter = opts.outIter;
outPmat = opts.outPmat;

% Create save covariance flag
if strcmp(outPmat,"all")
    savePs = true;
else
    savePs = false;
end 

% Create save output flag and preallocate
if strcmp(outIter,"all")
    saveAll = true;

    % Filtered Solution
    filSol.X = zerosCell([maxIter,1],[M,1]);
    filSol.dx = zerosCell([maxIter,1],[M,1]);
    filSol.bs = zerosCell([maxIter,1],[N,K]);

    % Smoothed Solution
    smoSol.X = zerosCell([maxIter,1],[M,1]);
    smoSol.dx = zerosCell([maxIter,1],[M,1]);

    % Covariance Matrices
    if savePs
        filSol.P = zerosCell([maxIter,1],[M,M]);
        smoSol.P = zerosCell([maxIter,1],[M,M]);
    end

    % Additional Outputs
    addOut.iter = zeros(maxIter,1);
    addOut.maxDiff = zeros(maxIter,1);
else
    % If not, no need, just use as a flag to save or not
    saveAll = false;
end 

% Check if timestep is constant
[cM,dt] = isuniform(t);

I_m = eye(M);

% If dt is not constant
if ~cM
    % Then Phi is nonconstant and we need to preallocate
    dt = diff(t);
end

Phis = zeros(M,M,N);
Phis(:,:,1) = I_m;

% Initialize filter with predicted state iteration and initial difference
X_pred0 = X_est0 - dx_est0;
iter = 0;
maxDiff = Inf;
numInc = 0;

while  maxDiff > tol && iter < maxIter && numInc < maxInc

    % Different calculations depending on how Q is calculated
    switch Qtype
        case "SNC"
            % TODO SNC Qtype

        case "DMC"
            % TODO DMC Qtype

        otherwise % Qs can be calculated individually

            % Set initial conditions for combined Phi and X_pred integration
            XpredPhi0 = [X_pred0;reshape(Phis(:,:,1),[],1)];

            % Update reference trajectory and STM
            [~,XpredsPhis] = vsIntFun(dotFunc,t,XpredPhi0,odeOpts);

            % Extract
            X_preds = XpredsPhis(:,1:M);
            rPhis = XpredsPhis(:,M+1:M+M*M);

            % Update to Incremental Phis
            for i = 2:N
                Phis(:,:,i) = reshape(rPhis(i,:),M,M)/reshape(rPhis(i-1,:),M,M);
            end

            if Qfunc
                % Evaluate Q matrix at each time as a function of time and state
                for i = 2:N
                    Qs(:,:,i) = Qdata(t(i),dt(i-1),X_preds(i,:));
                end
            end
    end

    if Rfunc
        % Evaluate R matrix at each time as function of time and state
        for i = 2:N
            Rs(:,:,i) = Rdata(t(i),dt(i-1),X_preds(i,:));
        end
    end

    if Hfunc
        % Evaluate H matrix at each point of the reference trajectory
        for i = 1:N
            Hs(:,:,i) = Hdata(t(i),X_preds(i,:));
        end
    end
    
    % Linearized Kalman Filter
    [X_ests,dx_ests,P_ests,P_preds,bs] = LKFpre(X_est0,dx_est0,P_0,Hs,Qs,Rs,zs,X_preds,Phis);

    % Rauch–Tung–Striebel Smoother
    [X_sms,dx_sms,P_sms] = RTSpre(X_preds,dx_ests,Phis,P_ests,P_preds);

    % Increment iterations
    iter = iter + 1;

    % Store Previous Difference
    pMaxDiff = maxDiff;

    % Calculate Difference
    maxDiff = max((abs(X_sms(1,:).'-X_pred0))./order);

    % Update numInc
    if maxDiff > pMaxDiff
        numInc = numInc + 1;
    else
        numInc = 0;
    end

    % Update
    X_pred0 = X_sms(1,:).';

    % Add the current iteration to the output
    if saveAll
        % Filtered Solution
        filSol.X{iter} = X_ests;
        filSol.dx{iter} = dx_ests;
        filSol.bs{iter} = bs;

        % Smoothed Solution
        smoSol.X{iter} = X_sms;
        smoSol.dx{iter} = dx_sms;

        % Covariance Matrices
        if savePs
            filSol.P{iter} = P_ests;
            smoSol.P{iter} = P_sms;
        end

        % Additional Outputs
        addOut.iter(iter) = iter;
        addOut.maxDiff(iter) = maxDiff;
    end
end

% Save only the last iteration 
if ~saveAll 
    % Filtered Solution
    filSol.X = X_ests;
    filSol.dx = dx_ests;
    filSol.bs = bs;

    % Smoothed Solution
    smoSol.X = X_sms;
    smoSol.dx = dx_sms;

    % Covariance Matrices
    if savePs
        filSol.P = P_ests;
        smoSol.P = P_sms;
    end

    % Additional Outputs
    addOut.iter = iter;
    addOut.maxDiff = maxDiff;
end

% Define Exit Flags and remove extra output
if maxDiff < tol
    addOut.exitState = "Filter Converged";
    if saveAll
        % Filtered Solution
        filSol.X(iter+1:end) = [];
        filSol.dx(iter+1:end) = [];
        filSol.bs(iter+1:end) = [];

        % Smoothed Solution
        smoSol.X(iter+1:end) = [];
        smoSol.dx(iter+1:end) = [];

        % Covariance Matrices
        if savePs
            filSol.P(iter+1:end) = [];
            smoSol.P(iter+1:end) = [];
        end

        % Additional Outputs
        addOut.iter(iter+1:end) = [];
        addOut.maxDiff(iter:end) = [];
    end
else
    if iter == maxIter
        addOut.exitState = "Maximum Iterations Reached";
    else
        addOut.exitState = "Filter Diverged";
    end
end

end % LKF_RTSpre Function

%% Helper Functions
function output = zerosCell(cellSize,matSize)
    %zeroCells Preallocates a cell array (k-rows, 1 column) with
    % INPUTS:
    % cellSize - [m1, m2, ... mk] Vector of cell array dimension lengths
    %  matSize - [n1, n2, ... nk] Vector of zero array dimension lengths
    % OUTPUT:
    % zeroCells - [m1, m2, ... mk] Cell array of [n1, n2, ... nk] zero matrices
    
    output = cell(cellSize);
    output(:) = {zeros(matSize)};
end

%% Validation Functions
function mustBe2Dor3DArrayorFunc(A)
% Checks if A is a 2D or 3D array or a function
if isa(A,"function_handle")
    if ~(isfile(func2str(A) + ".m") || exist(func2str(A),"builtin") == 5)
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

function mustBeValidFunction(A)
% Checks if A is a valid function handle
    if isa(A,"function_handle")
        if ~(isfile(func2str(A) + ".m") || ...
                exist(func2str(A),"file") == 2 || ...
                exist(func2str(A),"builtin") == 5)
            eid = "Type:notAValidFunctionHandle";
            msg = "Must be a valid function handle. Check for typos and " + ...
            "ensure that file is on your path.";
            error(eid,msg)
        end
    end
end

function observationsMustMatchTimeLength(t,obs)
% Checks if there are observations for each time step 
    if ~(size(obs,1) == length(t))
        eid = "Size:not2Dor3DArray";
        msg = "Must be a 2D or 3D array of doubles.";
        error(eid,msg)
    end
end

function mustBeValidOdeset(A)
% Checks if fields are consistent with odeset
    if ~isequal(fieldnames(odeset()),fieldnames(A))
        eid = "Fields:notMatchingOdeset";
        msg = "Ensure ode options are created using odeset before changing" ...
            + " values.";
        error(eid,msg)
    end
end

function mustBeValidInitConds(initConds)
% Checks if fields are consistent with initial conditions
% TODO add fieldname validation
    mustBeColumn(initConds.X_est0)
    mustBeColumn(initConds.dx_est0)
    mustBeMatrix(initConds.P_0)
    mustBeDouble(initConds.X_est0)
    mustBeDouble(initConds.dx_est0)
    mustBeDouble(initConds.P_0)
end

function mustBeValidDataMats(dataMats)
% Checks if fields are consistent with data matrices
% TODO add fieldname validation
    mustBe2Dor3DArrayorFunc(dataMats.Qdata)
    mustBe2Dor3DArrayorFunc(dataMats.Rdata)
    mustBe2Dor3DArrayorFunc(dataMats.Hdata)
    mustBeTextScalar(dataMats.Qtype)
    mustBeTextScalar(dataMats.Rtype)
    mustBeTextScalar(dataMats.Htype)
    mustBeMember(dataMats.Qtype,["const","data","SNC","DMC","func"])
    mustBeMember(dataMats.Rtype,["const","data","func"])
    mustBeMember(dataMats.Htype,["const","data","func"])
end

function mustBeValidSysFuncs(sysFuncs)
% Checks if fields are consistent with system functions
% TODO add fieldname validation
    mustBeValidFunction(sysFuncs.XdotPhidot)
    if nargin(sysFuncs.XdotPhidot) == 3
        if ~isfield(sysFuncs,"params")
            eid = "Fields:missingParametersField";
            msg = "State function needs parameters specified.";
            error(eid,msg)
        end
    else
        mustBeFunctionofTimeAndState(sysFuncs.XdotPhidot)
        if ~isfield(sysFuncs,"params")
            wid = "Fields:specifiedParametersUnused";
            msg = "Specified state function does not require parameters.";
            warning(wid,msg)
        end
    end
    mustBeValidFunction(sysFuncs.vsIntFun)
    mustBeOdeFunction(sysFuncs.vsIntFun)
    mustBeValidOdeset(sysFuncs.odeOpts)
    if isfield(sysFuncs,"params")
        if ~isstruct(sysFuncs.params)
            eid = "Type:paramsMustBeStruct";
            msg = "System function parameters must be specified in a " + ...
                "structure.";
            error(eid,msg)
        end
        if nargin(sysFuncs.XdotPhidot) == 2
            if ~isfield(sysFuncs,"params")
                wid = "Fields:specifiedParametersUnused";
                msg = "Specified parameters unused by state function.";
                warning(wid,msg)
            end           
        elseif nargin(sysFuncs.XdotPhidot) ~= 3
            eid = "Function:mustHaveThreeInputs";
            msg = "System function must accept three input arguments:" + ... 
            " time, state, and parameters.";
            error(eid,msg)
        end     
    end
end

function mustBeValidOpts(opts)
% Checks if fields are consistent with solving options
% TODO add fieldname validation
    mustBeDouble(opts.tol)
    mustBeNumeric(opts.maxIter)
    mustBeNumeric(opts.maxInc)
    mustBeScalar(opts.tol)
    mustBeScalar(opts.maxIter)
    mustBeScalar(opts.maxInc)
    mustBeFinite(opts.tol)
    mustBeFinite(opts.maxIter)
    mustBeFinite(opts.maxInc)
    mustBeNonempty(opts.tol)
    mustBeNonempty(opts.maxIter)
    mustBeNonempty(opts.maxInc)
    mustBeTextScalar(opts.outIter)
    mustBeTextScalar(opts.outPmat)
    mustBeMember(opts.outIter,["all","last"])
    mustBeMember(opts.outPmat,["all","none"])
end

function mustBeDouble(A)
% Checks if input is a double
    if ~isa(A,"double")
        eid = "Type:mustBeDouble";
        msg = "Variable must be a double.";
        error(eid,msg)
    end
end

function mustBeScalar(A)
    % Checks if input is a double
    if ~isscalar(A)
        eid = "Size:mustBeScalar";
        msg = "Variable must be a scalar.";
        error(eid,msg)
    end
end

function mustBeFunctionofTimeAndState(A)
    %Checks that the function accepts two input arguments
    if nargin(func2str(A)) ~= 2
        eid = "Function:mustHaveTwoInputs";
        msg = "State Functions must accept only time and state as " + ...
            "required arguments.";
        error(eid,msg)
    end
end

function mustBeOdeFunction(A)
    %Checks that the function matches MATLAB ode suite input arguments
    if nargin(func2str(A)) ~= -5
        eid = "Function:mustHaveMatlabOdeStyleInputs";
        msg = "Variable step integrator must have MATLAB style inputs.";
        error(eid,msg)
    end
end