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
geneFig = false;    % Whether or not to generate figures
saveFig = false;    % Whether or not to save figures
saveRes = false;    % Whether or not to save results
runCalc = true;    % Whether or not to run calculation
closeAfterSave = false; % Whether or not to close figures after saving
figfmt = "png"; % Format to save figures
% Recommended: "fig" or "mfig", "jpeg" or "png", "pdf" or "svg"

% RNG Setting for reproducibility
seed = 0;              % Random number seed
generator = "twister"; % Random number generator
rng(seed,generator)

if runCalc 
    % Parameters of Earth
    muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]
    rE = 6378.1363;    % Earth radius [km]
    wE = (2*pi)/86164.0905308; % Earth angular velocity (about z) [rad/s]

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
    XdotNL = @stateEQ_drag_9s6p;
    XdotPhidot = @stateSTM_drag_9s6p;
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
    sigmaMod = [1e-2;1e-1;1e-1]; % Standard Deviation of Initial Model Estimate
    sigmaPos = 1.0e-1*ones(3,1); % Standard Deviation of Initial Position Estimate
    sigmaVel = 1.0e-2*ones(3,1); % Standard Deviation of Initial Velocity Estimate
    covMod = sigmaMod.^2;        % Variance of Model Estimate
    covPos = sigmaPos.^2;        % Variance of Position Estimate
    covVel = sigmaVel.^2;        % Variance of Velocity Estimate

    % For H matrix
    states = ["rho_0","h_0","H","r_x","r_y","r_z","v_x","v_y","v_z"];
    nStates = length(states);
    observedStates = ["r_x","r_y","r_z"];
    obs = length(observedStates);

    % Set Matrix Types
    Qtype = "const"; % const, data, SNC, DMC, func
    Rtype = "const"; % const, data, func
    Htype = "const"; % const, data, func

    % Set Data Matrices
    Qdata = 1e-6*[covMod;covPos;covVel].*speye(nStates);
    Rdata = covMeas.*speye(obs);
    Hdata = [zeros(3),eye(3),zeros(3)];
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
    [satScen,initElems,initStates,initParams] = initialOrbits(constellation);

    % Use Initial Orbits to Prepare for Dataset Generation
    tf = numPeriods*initElems(1,7); % Final Time
    t = (t0:dt:tf)';                % Time Vector
    Nstep = length(t);              % Number of time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Initial Trajectory and Observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set all Parameters
    params.muE = muE;
    params.rE = rE;
    params.wE = wE;
    params.Cd = Cd;
    params.A = A;
    params.m = m;

    % Set constant parameters
    tparams.muE = params.muE;
    tparams.rE = params.rE;
    tparams.wE = params.wE;

    % Preallocate Nonlinear Trajectories
    Xnl = zeros(Nstep,nStates,nSats);

    % Calculate Nonlinear Trajectories
    for i = nSats:-1:1
        X0 = initStates(i,:).';
        tparams.Cd = params.Cd(i);
        tparams.A = params.A(i);
        tparams.m = params.m(i);
        [~,Xnl(:,:,i)] = vsIntFun(@(t,X) XdotNL(t,X,tparams),t,X0,odeOpts);
    end

    % Create Observations
    zs = Xnl(:,4:6,:) + sigmaMeas.'.*randn(Nstep,obs,nSats);

    % Set Initial Conditions
    X_est0 = squeeze(Xnl(1,:,:)).' + [sigmaMod;sigmaPos;sigmaVel].'.*randn(nSats,nStates);
    dx_est0 = zeros(nSats,nStates);
    P_0 = diag([covMod;covPos;covVel]).*ones(1,1,nSats);
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

    t1 = tic;
    % Perform Filtering for each satellite
    parfor i = 1:nSats % can change to parfor for now
        % In the future, we will need to assemble the matrices and state vector
        fprintf("Calculating LKF_RTSpre Pass for Satellite %d ...\n",nSats-i+1)
        % Set Initial Conditions
        initConds = struct();
        initConds.X_est0 = X_est0(i,:).';
        initConds.dx_est0 = dx_est0(i,:).';
        initConds.P_0 = P_0(:,:,i);
        t2 = tic;
        [filtered{i},smoothed{i},add_info{i}] = ...
            LKF_RTSpre(t,zs(:,:,i),initConds,dataMats,sysFuncs,opts);
        cTime = toc(t2);
        fprintf("Finished LKF_RTSpre Pass for Satellite %d in %.3f seconds\n",...
            nSats-i+1,cTime)
    end
    satSol = struct("filtered",filtered,"smoothed",smoothed,"add_info",add_info);
    tTime = toc(t1);
    clc;
    fprintf("Finished LKF_RTSpre Passes for all Satellites in %02d:%02d:%06.3f\n\n",...
        floor(tTime/3600),floor(mod(tTime,3600)/60),mod(tTime,60))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare Folder to Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveString = sprintf("%s_%s-%d_%s-%s-%s_%sIter_%sPmat_%.2EPer_dt%.1Es_Q-"+ ...
        "%s_R-%s_H-%s.mat",satNames,func2str(satFun),nSats,func2str(vsIntFun), ...
        func2str(XdotNL),func2str(XdotPhidot),outIter,outPmat,numPeriods,dt,...
        Qtype,Rtype,Htype);

    if ~isfolder("results")
        mkdir("results")
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear Redundant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear satFun inclin nPlanes phasing argLat satNames hSats nSats
    clear muE rE wE A Cd m rhoFun
    clear Qdata Rdata Hdata Qtype Rtype Htype QRH
    clear X0 XdotPhidot vsIntFun odeOpts
    clear smoTol maxIter maxInc outIter outPmat
    clear temp t1 t2 cTime BigSave resultString
    clear filtered smoothed add_info % These are the biggest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveRes && ~isfile("results/" + saveString)
    fprintf("Saving Workspace to file...\n")
    tic;
    save("results/"+saveString,"-v7.3")
    cTime = toc;
    fprintf("Results Saved in %02d:%02d:%06.3f\n\n", ...
        floor(cTime/3600),floor(mod(cTime,3600)/60),mod(cTime,60))
elseif saveRes && isfile("results/" + saveString)
    fprintf("File already exists, skipping saving.\n")
    fprintf("Delete or rename file below to save results:\n%s\n\n", ...
        saveString)
else
    fprintf("Result saving skipped, set ""saveRes = true""" + ...
        " to save calculation results.\n\n")
end
else
if saveRes
    fprintf("Saving results set to %s, but ...\n",string(saveRes))
    fprintf("Result calculation skipped, set ""runCalc = true""" + ...
        " to perform calculations.\n\n")
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~runCalc
    [fName,fLoc] = uigetfile({"*.mat","MAT-files"}, ...
        "Select a File","MultiSelect","off");
    if (ischar(fLoc) || isstring(fLoc)) && fLoc == fullfile(pwd,"/results/")
        if ischar(fLoc)
            fLoc = string(fLoc);
        end
        if ischar(fName)
            fName = string(fName);
        end
        load(fLoc+fName)
        fprintf("Loaded: %s\nFrom: %s\n",fName,fLoc)
        fileLoaded = true;
    else
        fprintf("Skipped loading. Load results from results subdirectory only.\n")
        fileLoaded = false;
    end
else
    fileLoaded = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if geneFig && (runCalc || fileLoaded)
    fprintf("\nBeginning Figure generation...\n")
    tic;
    nSats = constellation.nSats;
    rE = params.rE;
    if ~isfolder("images")
        mkdir("images")
    end

    % Smoothed trajectories
    close all;
    
    if saveFig
        st = figure;
    else
        figure;
    end
    [Xs,Ys,Zs] = sphere(20);
    surf(rE*Xs,rE*Ys,rE*Zs,"FaceColor",[0.3 0.3 0.3])
    hold on
    for i = 1:nSats
        if strcmp(opts.outIter,"last")
            %plot3(Xnl(:,4,i),Xnl(:,5,i),Xnl(:,6,i),":")
            plot3(satSol(i).smoothed.X(:,4), ...
                  satSol(i).smoothed.X(:,5), ...
                  satSol(i).smoothed.X(:,6))
        else
            %plot3(Xnl(:,4,i),Xnl(:,5,i),Xnl(:,6,i),":")
            plot3(satSol(i).smoothed.X{1}(:,4), ...
                  satSol(i).smoothed.X{1}(:,5), ...
                  satSol(i).smoothed.X{1}(:,6))
        end
    end
    hold off
    daspect([1 1 1])

    if saveFig
        saveas(st,"images/all_trajectories",figfmt)
    end

    % View Satellite Scenario
    ssv = satelliteScenarioViewer(satScen);

    % Constellation Visualization
    for sN = 1
        figure;
        surf(rE*Xs,rE*Ys,rE*Zs, ...
            "FaceColor",[0.3 0.3 0.3],"DisplayName","Earth")
        hold on
        plot3(Xnl(:,4,sN),Xnl(:,5,sN),Xnl(:,6,sN), ...
            "-k","DisplayName","Nonlinear")
        plot3(zs(:,4,sN),zs(:,5,sN),zs(:,6,sN), ...
            ".","DisplayName","Observations")
        if strcmp(opts.outIter,"last")
            plot3(satSol(sN).filtered.X(:,4), ...
                  satSol(sN).filtered.X(:,5), ...
                  satSol(sN).filtered.X(:,6), ...
                  "--","DisplayName","First Filter")
            plot3(satSol(sN).smoothed.X(:,4), ...
                  satSol(sN).smoothed.X(:,5), ...
                  satSol(sN).smoothed.X(:,6), ...
                  "--","DisplayName","First Smoother")
            plot3(satSol(sN).filtered.X(:,4), ...
                  satSol(sN).filtered.X(:,5), ...
                  satSol(sN).filtered.X(:,6), ...
                  ":","DisplayName","Last Filter")
            plot3(satSol(sN).smoothed.X(:,4), ...
                  satSol(sN).smoothed.X(:,5), ...
                  satSol(sN).smoothed.X(:,6), ...
                  ":","DisplayName","Last Smoother")
        else
            plot3(satSol(sN).filtered.X{1}(:,4), ...
                  satSol(sN).filtered.X{1}(:,5), ...
                  satSol(sN).filtered.X{1}(:,6), ...
                  "--","DisplayName","First Filter")
            plot3(satSol(sN).smoothed.X{1}(:,4), ...
                  satSol(sN).smoothed.X{1}(:,5), ...
                  satSol(sN).smoothed.X{1}(:,6), ...
                  "--","DisplayName","First Smoother")
            plot3(satSol(sN).filtered.X{end}(:,4), ...
                  satSol(sN).filtered.X{end}(:,5), ...
                  satSol(sN).filtered.X{end}(:,6),":", ...
                  "DisplayName","Last Filter")
            plot3(satSol(sN).smoothed.X{end}(:,4), ...
                  satSol(sN).smoothed.X{end}(:,5), ...
                  satSol(sN).smoothed.X{end}(:,6), ...
                  ":","DisplayName","Last Smoother")
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
        for i = (4:obs+3)
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
            saveas(fs(sN),sprintf("images/components%02d",sN),figfmt)
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
        for i = (4:obs+3)
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
            saveas(fd(sN),sprintf("images/error%02d",sN),figfmt)
        end
    end

    cTime = toc;
    fprintf("Finished Figure generation in %02d:%02d:%06.3f\n\n", ...
            floor(cTime/3600),floor(mod(cTime,3600)/60),mod(cTime,60))

    if saveFig && closeAfterSave
        close all
        delete(ssv);
        clear ssv
    end

    if (closeAfterSave && ~saveFig) && (isstring(figfmt) || ischar(figfmt))
        fprintf("Figures closing after saving set to %s, but...\n", ...
            string(closeAfterSave))
        fprintf("Figure saving as .%s is set to %s. ",figfmt,string(saveFig))
        fprintf("Set ""saveFig = true""" + ...
            " to save figures.\n")
    end
else
    if (runCalc || fileLoaded) && (isstring(figfmt) || ischar(figfmt))
        if closeAfterSave && saveFig
            fprintf("Figures saving as .%s set to %s ",figfmt,string(saveFig))
            fprintf("and closing after saving set to %s, but...\n", ...
                string(closeAfterSave))
        elseif closeAfterSave
            fprintf("Figures closing after saving set to %s, but...\n", ...
                string(closeAfterSave))
            fprintf("Figure saving as .%s set to %s. Additionally ...\n", ...
                figfmt,string(saveFig))
        elseif saveFig
            fprintf("Figure saving as .%s set to %s, but...\n", ...
                figfmt,string(saveFig))
        end
        fprintf("Figure generation skipped. Set ""geneFig = true""" + ...
                " to generate.\n")
    else
        fprintf("No calculations run or results loaded to plot.\n")
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear subdirectories to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clearPath
    rmpath("other_functions")
    rmpath("filter_smoother")
    rmpath("state_equations")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear Redundant Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ax fd fs h i lg st sN tiles tTime Xs Ys Zs rE nSats
clear saveString subDir nSats satData cTime figfmt resultString
clear closeAfterSave saveFig saveRes clearPath runCalc geneFig fileLoaded