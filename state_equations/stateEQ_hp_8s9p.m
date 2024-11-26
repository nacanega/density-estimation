function xyzdot = stateEQ_hp_8s9p(t,xyz,params)
%stateEQ_hp_8s6p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%      t - Current time [s]
%    xyz - [8,1] Current State
% params - Structure of parameters
%        .    muE - Earth gravitational parameter [km^3/s^2]
%        .     rE - Earth radius [km]
%        .     wE - Earth angular velocity [rad/s]
%        .     Cd - Drag Coefficient
%        .      A - Area [km^2]
%        .      m - Mass [kg]
%        .      n - Latitudinal Density Variation Exponent
%        .  UbFun - Diurnal Bulge Unit Vector
%        .     h0 - Base altitude [km]
% OUTPUT:
% xyzdot - [8,1] Current State Derivatives 

% Load parameters
muE = params.muE;
%rE = params.rE;
wE = params.wE;
Cd = params.Cd;
A = params.A;
m = params.m;
rhoFun = params.rhoFun; % Inputs: (t,[rs;vs;rho_min;rho_max],params)
%n = params.n;
%UbFun = params.UbFun; % Input: (t,[rs;vs],params); Output: [Ubx;Uby;Ubz]

% Angular velocity vector
wA = [0;0;wE];

% Preallocate
xyzdot = zeros(size(xyz));

% Positons
rs = xyz(1:3);
r = sqrt(rs.'*rs);

% Velocities
vs = xyz(4:6);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = sqrt(vRs.'*vRs);

% Density parameters
rho_min = xyz(7); 
rho_max = xyz(8);

% Densities
rho = rhoFun(t,[rs;vs;rho_min;rho_max],params);
%drhodh = -rho/H;

% State Derivatives
xyzdot(1:3) = vs;
xyzdot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;
%xyzdot(7:8) = zeros(2,1);

end