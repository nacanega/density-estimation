% Nolan Canegallo
% SingleSatSingleParam.m
% MAE 586 - Atmospheric Density Estimation Project
if exist("ssv","var") == 1
    delete(ssv)
end
clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subdirectories to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("state_equations")
addpath("state_equations/accelerations")
addpath("filter_smoother")
addpath("adensity_models")
addpath("other_functions")
% addpath("old")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of Earth
muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]
rE = 6378.1363;    % Earth radius [km]
wE = (2*pi)/86164.0905308; % Earth angular velocity (about z) [rad/s]
fE = 1/298.257223563;      % Earth flattening

% Problem Definition
nSats = 1;
altitudes = 781; % [km]
satType = "singleCircular";
satNames = "TestSat";
models.gravity = "";
models.density = "exp";
sysParams = ["C_D","A","m","h_0","H"];
satStates = ["r_x","r_y","r_z","v_x","v_y","v_z"];
modStates = "rho_0";
obsStates = ["r_x","r_y","r_z"];
combine = false;
nSatSt = length(satStates);
nModSt = length(modStates);
nStates = max(combine*nSats,1)*nSatSt + nModSt;
obs = length(obsStates);

% Parameters for drag (all satellites the same)
C_D = 2.2;
A = 1e1*1e-6;          % Cross sectional area [km^2]
m = 4;                 % Mass [kg]
h_0 = 700;             % Base altitude
H = 88.667;
defaultParams = [C_D A m h_0 H]; % TODO reimplement variations of these between sats

% Standard Deviations
sigmaR = 0.4e-3 * ones(3,1); % r_x, r_y, r_z
sigmaV = 1e-3 * ones(3,1); % v_x, v_y, v_z
sigmaS = [];%1e-2; % C_D
sigmaState = [sigmaR;sigmaV;sigmaS];
sigmaM = 1e-2;%;1e+1]; % rho_0, H
sigmaQ = [1e-6;1e-6;1e-6];
sigmaMeas = sigmaR; % r_x, r_y, r_z

% Covariances
P0 = 1e7;
covR = sigmaR.^2;
covV = sigmaV.^2;
covS = sigmaS.^2; % C_D
covM = sigmaM.^2;
covQ = sigmaQ.^2;
covState = sigmaState.^2;
covMeas = sigmaMeas.^2; 

% For Integration
t0 = 0;           % Initial Time [s]
numPeriods = 1; % Number of periods to calculate
dt = 1;           % Time Step [s]

% Set Matrix Types
Qtype = "SNC"; % const, data, SNC, DMC, func
Rtype = "const"; % const, data, func
Htype = "const"; % const, data, func

% Set Data Matrices
Pmati = P0*eye(nSatSt);
Pmatj = P0*eye(nModSt);
Rmati = covMeas.*eye(obs);
Hmati = [eye(obs),zeros(obs,nSatSt-obs)];
modMat = zeros(nModSt);
if strcmp(Qtype,"SNC")
    Qmati = covQ.*eye(obs);
    Qmatj = [];
    Qrotm = true; % true, false
else
    Qmati = covState.*eye(nSatSt);
    Qmatj = 1e-14*eye(nModSt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flags and Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*UNRCH> % Flag to suppress unreachable code warning
% Flags
clearPath = false; % Whether or not to clear path after running
geneFig = false;    % Whether or not to generate figures
saveFig = false;    % Whether or not to save figures
saveRes = false;    % Whether or not to save results
runCalc = true;    % Whether or not to run calculation
closeAfterSave = false; % Whether or not to close figures after saving
figfmt = "png"; % Format to save figures
% Recommended: "fig" or "mfig", "jpeg" or "png", "pdf" or "svg"

% Integrator
vsIntFun = @ode113;
odeOpts = odeset("AbsTol",1e-9,"RelTol",1e-8);

% RNG Setting for reproducibility
seed = 0;              % Random number seed
generator = "twister"; % Random number generator
rng(seed,generator)

% Set Filter Options
smoTol = 2e-3;
maxIter = 100;
maxInc = 10;
outIter = "all"; % all or last % WARNING all uses a lot of RAM
outPmat = "all"; % all or none

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Satellite Initial Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Initial Satellite Constellation Parameters
constellation = initialConstell(satType);

% Change options from the defaults
constellation.nSats = nSats;
constellation.semimajor = altitudes + rE;
constellation.satNames = satNames;

% Generate Satellite Initial States
[satScen,initElems,initStates,varParams] = ...
    initialStates( ...
        constellation,satStates,modStates,sysParams,defaultParams,true);

% Use Initial Orbits to Prepare for Dataset Generation
tf = numPeriods*initElems(end);   % Final Time
t = (t0:dt:tf)';                  % Time Vector
Nstep = length(t);                % Number of time steps
%odeOpts.MaxStep = maxStepPercent*(tf-t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Initial Trajectory and Observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Functions
[gravAccel,dragAccel,GammaROT,stateEQ,stateSTM] = ...
    autoEquations(satStates,modStates,models,nStates);

% Set parameters
params.muE = muE;
params.rE = rE;
params.wE = wE;
%params.fE = fE;
params.C_D = C_D;
params.A = A;
params.m = m;
%params.rho_0 = rho_0;
params.h_0 = h_0;
params.H = H;
params.gravFunc = gravAccel;
params.dragFunc = dragAccel;
params.nSatStates = nSatSt;
params.nModStates = nModSt;
params.nSats = nSats;

X0 = initStates;
orders = orderFun(X0,nSats,nSatSt,nModSt);

[~,Xnl] = vsIntFun(@(t,X) stateEQ(t,X,params),t,X0,odeOpts);

figure
[Xs,Ys,Zs] = sphere(20);
surf(rE*Xs,rE*Ys,rE*Zs,"FaceColor",[0.3 0.3 0.3])
hold on
plot3(Xnl(:,1),Xnl(:,2),Xnl(:,3))
hold off
daspect([1 1 1])

figure
tiledlayout(4,1)
nexttile
plot(Xnl(:,4))
nexttile
plot(Xnl(:,5))
nexttile
plot(Xnl(:,6))
nexttile
plot(Xnl(:,7))

zs = measNoise(Xnl,sigmaMeas,nSats,Hmati);
X_est0 = stateNoise(Xnl(1,:)',sigmaState,sigmaM,nSats,nSatSt);
dx_est0 = zeros(nStates,1);
P_0 = blockMat(nSats,Pmati,Pmatj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Data Matrices
if strcmp(Qtype,"SNC")
    dataMats.Qdata = blockMat(nSats,Qmati);
else
    dataMats.Qdata = blockMat(nSats,Qmati,Qmatj);
end
dataMats.Rdata = blockMat(nSats,Rmati);
dataMats.Hdata = blockMat(nSats,Hmati,modMat);
dataMats.Hdata(end-(nModSt-1):end,:) = [];
dataMats.Qtype = Qtype;
dataMats.Rtype = Rtype;
dataMats.Htype = Htype;
if strcmp(Qtype,"SNC")
    dataMats.Qrotm = Qrotm;
end

% Set System Functions
sysFuncs.XdotPhidot = stateSTM;
sysFuncs.vsIntFun = vsIntFun;
sysFuncs.odeOpts = odeOpts;
sysFuncs.gamRotFun = GammaROT;
sysFuncs.params = params;

% Set Remaining Filter Options
opts.tol = smoTol;
opts.maxIter = maxIter;
opts.maxInc = maxInc;
opts.outIter = outIter;
opts.outPmat = outPmat;
opts.order = orders;

% Initial Conditions
initConds.X_est0 = X_est0;
initConds.dx_est0 = dx_est0;
initConds.P_0 = P_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[filtered,smoothed,add_info] = LKF_RTSpre(t,zs,initConds,dataMats,sysFuncs,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
figure
[Xs,Ys,Zs] = sphere(20);
surf(rE*Xs,rE*Ys,rE*Zs,"FaceColor",[0.3 0.3 0.3])
hold on
plot3(zs(:,1),zs(:,2),zs(:,3),".")
plot3(filtered(1).X(:,1),filtered(1).X(:,2),filtered(1).X(:,3),".")
plot3(smoothed(1).X(:,1),smoothed(1).X(:,2),smoothed(1).X(:,3),".")
plot3(zs(1,1),zs(1,2),zs(1,3),"o")
plot3(filtered(1).X(end,1),filtered(1).X(end,2),filtered(1).X(end,3),"or")
plot3(filtered(1).X(1,1),filtered(1).X(1,2),filtered(1).X(1,3),"og")
hold off
daspect([1 1 1])

figure
semilogy(Xnl(:,end))
hold on
semilogy(abs(filtered(1).X(:,end)),".")
semilogy(abs(smoothed(1).X(:,end)),".")
semilogy(abs(filtered(end).X(:,end)),".")
semilogy(abs(smoothed(end).X(:,end)),".")
hold off
title("rho_0")

figure
t1 = tiledlayout(3,1);
title(t1,"r Vectors")

nexttile
plot(zs(:,1)-Xnl(:,1),".")
hold on
plot(filtered(1).X(:,1)-Xnl(:,1),".")
plot(smoothed(1).X(:,1)-Xnl(:,1),".")
plot(filtered(end).X(:,1)-Xnl(:,1),".")
plot(smoothed(end).X(:,1)-Xnl(:,1),".")
hold off

nexttile
plot(zs(:,2)-Xnl(:,2),".")
hold on
plot(filtered(1).X(:,2)-Xnl(:,2),".")
plot(smoothed(1).X(:,2)-Xnl(:,2),".")
plot(filtered(end).X(:,2)-Xnl(:,2),".")
plot(smoothed(end).X(:,2)-Xnl(:,2),".")
hold off

nexttile
plot(zs(:,3)-Xnl(:,3),".")
hold on
plot(filtered(1).X(:,3)-Xnl(:,3),".")
plot(smoothed(1).X(:,3)-Xnl(:,3),".")
plot(filtered(end).X(:,3)-Xnl(:,3),".")
plot(smoothed(end).X(:,3)-Xnl(:,3),".")
hold off

figure
t2 = tiledlayout(3,1);
title(t2,"v Vectors");

nexttile
hold on
plot(filtered(1).X(:,4)-Xnl(:,4),".")
plot(smoothed(1).X(:,4)-Xnl(:,4),".")
plot(filtered(end).X(:,4)-Xnl(:,4),".")
plot(smoothed(end).X(:,4)-Xnl(:,4),".")
hold off

nexttile

hold on
plot(filtered(1).X(:,5)-Xnl(:,5),".")
plot(smoothed(1).X(:,5)-Xnl(:,5),".")
plot(filtered(end).X(:,5)-Xnl(:,5),".")
plot(smoothed(end).X(:,5)-Xnl(:,5),".")
hold off

nexttile
hold on
plot(filtered(1).X(:,6)-Xnl(:,6),".")
plot(smoothed(1).X(:,6)-Xnl(:,6),".")
plot(filtered(end).X(:,6)-Xnl(:,6),".")
plot(smoothed(end).X(:,6)-Xnl(:,6),".")
hold off