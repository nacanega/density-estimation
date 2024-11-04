% Nolan Canegallo
% SimpleDensityEstimation.m
% MAE 586 - Atmospheric Density Estimation Project
clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subdirectories to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("state_equations")
addpath("filter_smoother")
addpath("other_functions")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*UNRCH> % Flag to suppress unreachable code warning
% Flags
clearPath = false; % Whether or not to clear path after running
saveFig = true; % Whether or not to save figures
saveRes = true; % Whether or not to save results
closeAfterSave = true; % Whether or not to close figures after saving
figfmt = "png"; % Format to save figures:
% Recommended: "fig" or "mfig", "jpeg" or "png", "pdf" or "svg"

% RNG Setting for reproducibility
seed = 0;              % Random number seed
generator = "twister"; % Random number generator
rng(seed,generator)

% Parameters of Earth
muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]
rE = 6378.1363;    % Earth radius [km]
wE = (2*pi)/86164.0905308; % Earth angular velocity (about z) [rad/s]
rhoFun = @expDensity;      % Atmospheric Density function

% Parameters for drag (all satellites the same)
Cd = 2.17;                 % Nominal Drag coefficient
A = 1e1*1e-6;                % Cross sectional area [km^2]
m = 4;                     % Mass [kg]

% For Satellite Constellation
satFun = @walkerStar; % Function to generate constellation
hSats = 781;          % Satellite Altitude [km]
nSats = 66;           % Number of Satellites % Saving everything ~40 GB
inclin = 86.4;        % Inclination of orbits [deg]
nPlanes = 6;          % Number of planes
phasing = 2;          % Phasing
argLat = 0;           % Argument of latitude [deg]
satNames = "Iridium"; % Name String

% For Integration
t0 = 0;           % Initial Time [s]
numPeriods = 1; % Number of periods to calculate
dt = 1;           % Time Step [s]

% Set System Functions
XdotNL = @stateEQ_drag_6s7p;
XdotPhidot = @stateSTM_drag_6s7p;
vsIntFun = @ode89;
odeOpts = odeset("AbsTol",1e-9,"RelTol",1e-8);

% Set Filter Options
smoTol = 1e-6;
maxIter = 100;
maxInc = 10;
outIter = "all"; % all or last % WARNING all uses a lot of RAM
outPmat = "all"; % all or none

% Matrix Options
% For R matrix
sigmaMeas = 0.4e-3*ones(3,1); % Standard Deviation of Position Measurement
covMeas = sigmaMeas.^2;       % Variance of Position Measurement

% For P matrix
sigmaPos = 1.0e-1*ones(3,1); % Standard Deviation of Initial Position Estimate
sigmaVel = 1.0e-2*ones(3,1); % Standard Deviation of Initial Velocity Estimate
covPos = sigmaPos.^2;        % Variance of Position Estimate
covVel = sigmaVel.^2;        % Variance of Velocity Estimate

% For H matrix
states = ["r_x","r_y","r_z","v_x","v_y","v_z"];
nStates = length(states);
observedStates = ["r_x","r_y","r_z"];
obs = length(observedStates);

% Set Matrix Types
Qtype = "const"; % const, data, SNC, DMC, func
Rtype = "const"; % const, data, func
Htype = "const"; % const, data, func

% Set Data Matrices
Qdata = 1e-6*[covPos;covVel].*speye(nStates);
Rdata = covMeas.*speye(obs);
Hdata = speye(obs,nStates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Satellite Initial Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Initial Satellite Constellation Parameters
constellation = initialConstell("test");

% Set Constellation parameters
constellation.satFun = satFun;
constellation.radius = hSats + rE;
constellation.inclin = inclin;
constellation.nSats = nSats;
constellation.nPlanes = nPlanes;
constellation.phasing = phasing;
constellation.argLat = argLat;
constellation.satNames = satNames;

% Generate Satellite Initial States
[satScen,initElems,initStates] = initialOrbits(constellation);

% Use Initial Orbits to Prepare for Dataset Generation
tf = numPeriods*initElems(1,7); % Final Time
t = (t0:dt:tf)';                % Time Vector
Nstep = length(t);              % Number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Initial Trajectory and Observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set All Parameters
params.muE = muE;
params.rE = rE;
params.wE = wE;
params.rhoFun = rhoFun;
params.Cd = Cd;
params.A = A;
params.m = m;

% Preallocate Nonlinear Trajectories
Xnl = zeros(Nstep,nStates,nSats);

% Calculate Nonlinear Trajectories
for i = nSats:-1:1
    X0 = initStates(i,:).';
    [~,Xnl(:,:,i)] = vsIntFun(@(t,X) XdotNL(t,X,params),t,X0,odeOpts);
end

% Create Observations
zs = Xnl(:,1:3,:) + sigmaMeas.'.*randn(Nstep,obs,nSats);

% Set Initial Conditions
X_est0 = squeeze(Xnl(1,:,:)).' + [sigmaPos;sigmaVel].'.*randn(nSats,nStates);
dx_est0 = zeros(nSats,nStates);
P_0 = diag([covPos;covVel]).*ones(1,1,nSats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Matrix Types for Preallocation
QRH.Qdata = Qdata;
QRH.Rdata = Rdata;
QRH.Hdata = Hdata;
QRH.Qtype = Qtype;
QRH.Rtype = Rtype;
QRH.Htype = Htype;

% Generate Initial Filter Data Structures
[initConds,dataMats,sysFuncs] = initialStructs(obs,nStates,nSats,Nstep,QRH);

% Set Data Matrices
dataMats.Qdata = Qdata;
dataMats.Rdata = Rdata;
dataMats.Hdata = Hdata;
dataMats.Qtype = Qtype;
dataMats.Rtype = Rtype;
dataMats.Htype = Htype;

% Set System Functions
sysFuncs.XdotPhidot = XdotPhidot;
sysFuncs.vsIntFun = vsIntFun;
sysFuncs.odeOpts = odeOpts;
sysFuncs.params = params;

% Set Remaining Filter Options
opts.tol = smoTol;
opts.maxIter = maxIter;
opts.maxInc = maxInc;
opts.outIter = outIter;
opts.outPmat = outPmat;

% Setup structure to hold solution
filtered = cell(nSats,1);
smoothed = cell(nSats,1);
add_info = cell(nSats,1);

tic;
% Perform Filtering for each satellite
parfor i = 1:nSats % can change to parfor for now
    % In the future, we will need to assemble the matrices and state vector
    fprintf("Calculating LKF_RTSpre Pass for Satellite %d ...\n",nSats-i+1)
    % Set Initial Conditions
    initConds = struct();
    initConds.X_est0 = X_est0(i,:).';
    initConds.dx_est0 = dx_est0(i,:).';
    initConds.P_0 = P_0(:,:,i);
    tic;
    [filtered{i},smoothed{i},add_info{i}] = ...
        LKF_RTSpre(t,zs(:,:,i),initConds,dataMats,sysFuncs,opts);
    cTime = toc;
    fprintf("Finished LKF_RTSpre Pass for Satellite %d in %.3f seconds\n",...
        nSats-i+1,cTime)
end
satSol = struct("filtered",filtered,"smoothed",smoothed,"add_info",add_info);
tTime = toc;
clc;
fprintf("Finished LKF_RTSpre Passes for all Satellites in: %02d:%06.3f\n\n",...
    floor(tTime/60),mod(tTime,60))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveString = sprintf("%s_%s-%d_%s-%s-%s_%sIter_%sPmat_%.2EPer_dt%.1Es_Q-"+ ...
    "%s_R-%s_H-%s.mat",satNames,func2str(satFun),nSats,func2str(vsIntFun), ...
    func2str(XdotNL),func2str(XdotPhidot),outIter,outPmat,numPeriods,dt,...
    Qtype,Rtype,Htype);

if ~isfolder("results")
    mkdir("results")
end

subDir = sprintf("results\\%s",satNames);
if ~isfolder(subDir)
    mkdir(subDir)
end

if saveRes && ~isfile("results\\" + saveString)
    fprintf("Saving Workspace to file...\n")
    tic;
    save("results\\"+saveString,"-v7.3")
    for i = 1:nSats
        fprintf("Saving Satellite %d Results to file...\n",i)
        resultString = subDir + sprintf("\\satellite%d.mat",i);
        temp = satSol(i);
        save(resultString,"temp","-v7.3");
    end
    cTime = toc;
    fprintf("Results Saved in %02d:%06.3f\n\n",floor(cTime/60),mod(cTime,60))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear Redundant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear satFun inclin nSats nPlanes phasing argLat satNames
clear Qdata Rdata Hdata Qtype Rtype Htype QRH
clear X0 XdotPhidot vsIntFun odeOpts
clear smoTol maxIter maxInc outIter outPmat
clear temp cTime BigSave subDir saveString resultString
clear filtered smoothed add_info % These are the biggest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load("Iridium_walkerStar-66_ode89-stateEQ_nodrag-stateSTM_nodrag_allIter_"+ ...
%   "allPmat_1.00E+00Per_dt1.0E+00s_Q-const_R-const_H-const.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSats = constellation.nSats;
if ~isfolder("images")
    mkdir("images")
end

% Smoothed trajectories
close all;

st = figure;
[Xs,Ys,Zs] = sphere(20);
surf(rE*Xs,rE*Ys,rE*Zs,"FaceColor",[0.3 0.3 0.3])
hold on
for i = 1:nSats
    if strcmp(opts.outIter,"last")
        %plot3(Xnl(:,1,i),Xnl(:,2,i),Xnl(:,3,i),":")
        plot3(satSol(i).smoothed.X(:,1),satSol(i).smoothed.X(:,2),...
            satSol(i).smoothed.X(:,3))
    else
        %plot3(Xnl(:,1,i),Xnl(:,2,i),Xnl(:,3,i),":")
        plot3(satSol(i).smoothed.X{1}(:,1),satSol(i).smoothed.X{1}(:,2),...
            satSol(i).smoothed.X{1}(:,3))
    end
end
hold off
daspect([1 1 1])

if saveFig
    saveas(st,"images\\all_trajectories",figfmt)
end
% View Satellite Scenario

ssv = satelliteScenarioViewer(satScen);

% Constellation Visualization
for sN = 1
    figure;
    surf(rE*Xs,rE*Ys,rE*Zs,"FaceColor",[0.3 0.3 0.3],"DisplayName","Earth")
    hold on
    plot3(Xnl(:,1,sN),Xnl(:,2,sN),Xnl(:,3,sN),"-k","DisplayName","Nonlinear")
    plot3(zs(:,1,sN),zs(:,2,sN),zs(:,3,sN),".","DisplayName","Observations")
    if strcmp(opts.outIter,"last")
        plot3(satSol(sN).filtered.X(:,1),satSol(sN).filtered.X(:,2), ...
            satSol(sN).filtered.X(:,3),"--","DisplayName","First Filter")
        plot3(satSol(sN).smoothed.X(:,1),satSol(sN).smoothed.X(:,2), ...
            satSol(sN).smoothed.X(:,3),"--","DisplayName","First Smoother")
        plot3(satSol(sN).filtered.X(:,1),satSol(sN).filtered.X(:,2), ...
            satSol(sN).filtered.X(:,3),":","DisplayName","Last Filter")
        plot3(satSol(sN).smoothed.X(:,1),satSol(sN).smoothed.X(:,2), ...
            satSol(sN).smoothed.X(:,3),":","DisplayName","Last Smoother")
    else
        plot3(satSol(sN).filtered.X{1}(:,1),satSol(sN).filtered.X{1}(:,2), ...
            satSol(sN).filtered.X{1}(:,3),"--","DisplayName","First Filter")
        plot3(satSol(sN).smoothed.X{1}(:,1),satSol(sN).smoothed.X{1}(:,2), ...
            satSol(sN).smoothed.X{1}(:,3),"--","DisplayName","First Smoother")
        plot3(satSol(sN).filtered.X{end}(:,1),satSol(sN).filtered.X{end}(:,2), ...
            satSol(sN).filtered.X{end}(:,3),":","DisplayName","Last Filter")
        plot3(satSol(sN).smoothed.X{end}(:,1),satSol(sN).smoothed.X{end}(:,2), ...
            satSol(sN).smoothed.X{end}(:,3),":","DisplayName","Last Smoother")
    end
    hold off
    xlabel(observedStates(1)+", km")
    ylabel(observedStates(2)+", km")
    zlabel(observedStates(3)+", km")
    title(sprintf("Satellite %d",sN))
    legend show
    grid on
    daspect([1 1 1])
end

% Trajectory Visualizations
if saveFig && closeAfterSave
    close all;
end

fs = zeros(nSats,1);
for sN = 1:nSats
    fs(sN) = figure;
    h = gobjects(6,3);
    tiles = tiledlayout(3,1);
    title(tiles,sprintf("Satellite %d",sN))
    for i = 1:obs
        ax = nexttile;
        h(1,i) = plot(t,Xnl(:,i,sN),"-k");
        hold on
        h(2,i) = plot(t,zs(:,i,sN),".");
        if strcmp(opts.outIter,"last")
            h(3,i) = plot(t,satSol(sN).filtered.X(:,i),"--");
            h(4,i) = plot(t,satSol(sN).smoothed.X(:,i),"--");
            h(5,i) = plot(t,satSol(sN).filtered.X(:,i),":");
            h(6,i) = plot(t,satSol(sN).smoothed.X(:,i),":");
        else
            h(3,i) = plot(t,satSol(sN).filtered.X{1}(:,i),"--");
            h(4,i) = plot(t,satSol(sN).smoothed.X{1}(:,i),"--");
            h(5,i) = plot(t,satSol(sN).filtered.X{end}(:,i),":");
            h(6,i) = plot(t,satSol(sN).smoothed.X{end}(:,i),":");
        end
        hold off
        xlim([t(1) t(end)])
        xlabel("Time, s")
        ylabel(observedStates(i)+", km")
        grid on
    end
    lg = legend(ax,h(:,3),["Nonlinear Solution","Observations", ...
        "First Filter Iteration","First Smoother Iteration", ...
        "Last Filter Iteration","Last Smoother Iteration"]);
    lg.Layout.Tile = "East";
    if saveFig
        saveas(fs(sN),sprintf("images\\components%02d",sN),figfmt)
    end
end

if saveFig && closeAfterSave
    close all;
end

% Error Visualizations

fd = zeros(nSats,1);

for sN = 1:nSats
    fd(sN) = figure;
    h = gobjects(5,3);
    tiles = tiledlayout(3,1);
    title(tiles,sprintf("Satellite %d",sN))
    for i = 1:obs
        ax = nexttile;
        h(1,i) = plot(t,zs(:,i,sN)-Xnl(:,i,sN),".");
        hold on
        if strcmp(opts.outIter,"last")
            h(2,i) = plot(t,satSol(sN).filtered.X(:,i)-Xnl(:,i,sN),".");
            h(3,i) = plot(t,satSol(sN).smoothed.X(:,i)-Xnl(:,i,sN),".");
            h(4,i) = plot(t,satSol(sN).filtered.X(:,i)-Xnl(:,i,sN),".");
            h(5,i) = plot(t,satSol(sN).smoothed.X(:,i)-Xnl(:,i,sN),".");
        else
            h(2,i) = plot(t,satSol(sN).filtered.X{1}(:,i)-Xnl(:,i,sN),".");
            h(3,i) = plot(t,satSol(sN).smoothed.X{1}(:,i)-Xnl(:,i,sN),".");
            h(4,i) = plot(t,satSol(sN).filtered.X{end}(:,i)-Xnl(:,i,sN),".");
            h(5,i) = plot(t,satSol(sN).smoothed.X{end}(:,i)-Xnl(:,i,sN),".");
        end
        hold off
        xlim([t(1) t(end)])
        ylim(3*[-sigmaMeas(1) sigmaMeas(1)])
        xlabel("Time, s")
        ylabel(observedStates(i)+" Error, km")
        grid on
    end
    lg = legend(ax,h(:,3),["Observations", ...
        "First Filter Iteration","First Smoother Iteration", ...
        "Last Filter Iteration","Last Smoother Iteration"]);
    lg.Layout.Tile = "East";
    if saveFig
        saveas(fd(sN),sprintf("images\\error%02d",sN),figfmt)
    end
end

if saveFig && closeAfterSave
    close all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear Redundant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ax fd fs h i lg st sN tiles tTime Xs Ys Zs
clear muE rE wE
clear closeAfterSave saveFig saveRes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear subdirectories to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clearPath
    rmpath("other_functions")
    rmpath("filter_smoother")
    rmpath("state_equations")
end