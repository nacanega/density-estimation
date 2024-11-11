function xyzdot = stateEQ_phasedrag_9s8p(~,xyz,params)
%stateSTM_phasedrag_9s8p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyzPhi - [9,1] Current State and STM
%    params - Structure of parameters
%           .    muE - Earth gravitational parameter [km^3/s^2]
%           .     rE - Earth radius [km]
%           .     wE - Earth angular velocity [rad/s]
%           . rhoFun - Density function handle
%           .     Cd - Drag Coefficient
%           .      A - Area [km^2]
%           .      m - Mass [kg]
%           . lambda - Sun parameter [?]
%           
% OUTPUT:
% xyzPhidot - [9,1] Current State and STM Derivatives 

% Load parameters
muE = params.muE;
rE = params.rE;
wE = params.wE;
rhoFun = params.rhoFun;
Cd = params.Cd;
A = params.A;
m = params.m;
lambda = params.lambda;

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
vR = norm(vRs);

% Density
rho = rhoFun(~,r,lambda);

% State Derivatives
xyzdot(1:3) = vs;
xyzdot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;

end