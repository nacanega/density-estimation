% Nolan Canegallo
% SimpleTest.m
% MAE 586 - Atmospheric Density Estimation Project
% 28 October 2024
clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(13)

cov = 1e-2;
sigma = sqrt(cov);
g = 9.8;

t0 = 0;
tf = 3;
dt = 0.01;
t = (t0:dt:tf).';
N = length(t);

X0 = [0;1.5*g;g;g];
X_est0 = [X0(1:2);10;10] + [-1;1;randn(1);randn(1)];

XdotNL = @stateEQ_ffdrag;
XdotPhidot = @stateSTM_ffdrag;
vsIntFun = @ode45;
odeOpts = odeset();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear Solutions and Observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,Xnl] = vsIntFun(XdotNL,t,X0,odeOpts);
[~,Xrf] = vsIntFun(XdotNL,t,X_est0,odeOpts);
zs = Xnl(:,1:2) + sigma*randn(N,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
dataMats.Qdata = (1)^2*dt^2*[dt^2/4*eye(2),dt/2*eye(2);dt/2*eye(2),eye(2)];
dataMats.Rdata = cov*eye(2);
dataMats.Hdata = eye(2,4);
dataMats.Qtype = "const";
dataMats.Rtype = "const";
dataMats.Htype = "const";

% Set System Functions
sysFuncs.XdotPhidot = XdotPhidot;
sysFuncs.vsIntFun = vsIntFun;
sysFuncs.odeOpts = odeOpts;

% Set Remaining Filter Options
opts.tol = 1e-6;
opts.maxIter = 150;
opts.maxInc = 10;
opts.outIter = "all";
opts.outPmat = "all";

initConds.X_est0 = X_est0;
initConds.dx_est0 = zeros(4,1);
initConds.P_0 = 1/cov*eye(4);

[filtered,smoothed,add_info] = LKF_RTSpre(t,zs,initConds,dataMats,sysFuncs,opts);

% Second Reference Trajectory
Xrf20 = smoothed.X{1}(1,:);
[~,Xrf2] = vsIntFun(XdotNL,t,Xrf20,odeOpts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
TLR = tiledlayout(2,3);

nexttile([2 2])
plot(Xnl(:,1),Xnl(:,2),LineWidth=1.2)
hold on
plot(zs(:,1),zs(:,2),".")
plot(Xrf(:,1),Xrf(:,2),LineWidth=1.2)
plot(filtered.X{1}(:,1),filtered.X{1}(:,2),LineWidth=1.2)
plot(smoothed.X{1}(:,1),smoothed.X{1}(:,2),LineWidth=1.2)
plot(Xrf2(:,1),Xrf2(:,2),":",LineWidth=1.2)
plot(filtered.X{end}(:,1),filtered.X{end}(:,2),"--",LineWidth=1.2)
plot(smoothed.X{end}(:,1),smoothed.X{end}(:,2),"--",LineWidth=1.2)
hold off
grid on

legend(["NL","OBS","REF","LKF1","RTS1","REF2","LKFN","RTSN"], ...
    "Location","southoutside","Orientation","horizontal","NumColumns",4)
xlabel("X Position, m")
ylabel("Y Position, m")
title("Trajectories")

% X
nexttile
semilogy(t,abs(zs(:,1)-Xnl(:,1)),".")
hold on
semilogy(t,abs(Xrf(:,1)-Xnl(:,1)),".")
semilogy(t,abs(filtered.X{1}(:,1)-Xnl(:,1)),".")
semilogy(t,abs(smoothed.X{1}(:,1)-Xnl(:,1)),".")
semilogy(t,abs(filtered.X{end}(:,1)-Xnl(:,1)),".")
semilogy(t,abs(smoothed.X{end}(:,1)-Xnl(:,1)),".")
hold off
grid on

legend(["OBS","REF","LKF1","RTS1","LKFN","RTSN"],"Location","eastoutside")
ylabel("X Position Error, m")
title("Absolute Error in X Position")

% Y
nexttile
semilogy(t,abs(zs(:,2)-Xnl(:,2)),".")
hold on
semilogy(t,abs(Xrf(:,2)-Xnl(:,2)),".")
semilogy(t,abs(filtered.X{1}(:,2)-Xnl(:,2)),".")
semilogy(t,abs(smoothed.X{1}(:,2)-Xnl(:,2)),".")
semilogy(t,abs(filtered.X{end}(:,2)-Xnl(:,2)),".")
semilogy(t,abs(smoothed.X{end}(:,2)-Xnl(:,2)),".")
hold off
grid on

legend(["OBS","REF","LKF1","RTS1","LKFN","RTSN"],"Location","eastoutside")
xlabel("Time, s")
ylabel("Y Position Error, m")
title("Absolute Error in Y Position")

%% Velocity
figure
TLV = tiledlayout(2,2);

% VX
nexttile
plot(t,Xnl(:,3))
hold on
plot(t,Xrf(:,3),".")
plot(t,filtered.X{1}(:,3),".")
plot(t,smoothed.X{1}(:,3),".")
hold off
grid on

title("X Velocity")


% VY
nexttile
plot(t,Xnl(:,4))
hold on
plot(t,Xrf(:,4),".")
plot(t,filtered.X{1}(:,4),".")
plot(t,smoothed.X{1}(:,4),".")
hold off
grid on

legend(["NL","REF","LKF1","RTS1"],"Location","eastoutside")
title("Y Velocity")

% VX Error
nexttile
semilogy(t,abs(Xrf(:,3)-Xnl(:,3)),".")
hold on
semilogy(t,abs(filtered.X{1}(:,3)-Xnl(:,3)),".")
semilogy(t,abs(smoothed.X{1}(:,3)-Xnl(:,3)),".")
semilogy(t,abs(filtered.X{end}(:,3)-Xnl(:,3)),".")
semilogy(t,abs(smoothed.X{end}(:,3)-Xnl(:,3)),".")
hold off
grid on

title("X Velocity Absolute Error")


% VY Error
nexttile
semilogy(t,abs(Xrf(:,4)-Xnl(:,4)),".")
hold on
semilogy(t,abs(filtered.X{1}(:,4)-Xnl(:,4)),".")
semilogy(t,abs(smoothed.X{1}(:,4)-Xnl(:,4)),".")
semilogy(t,abs(filtered.X{end}(:,4)-Xnl(:,4)),".")
semilogy(t,abs(smoothed.X{end}(:,4)-Xnl(:,4)),".")
hold off
grid on

legend(["REF","LKF1","RTS1","LKFN","RTSN"],"Location","eastoutside")
title("Y Velocity Absolute Error")

% Shared Labels
xlabel(TLV,"Time, s")
ylabel(TLV,"Velocity, m/s")