function xyzdot = stateEQ_drag_9s6p(~,xyz,params)
%stateEQ_drag_9s6p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyz - [9,1] Current State
% params - Structure of parameters
%        .    muE - Earth gravitational parameter [km^3/s^2]
%        .     rE - Earth radius [km]
%        .     wE - Earth angular velocity [rad/s]
%        .     Cd - Drag Coefficient
%        .      A - Area [km^2]
%        .      m - Mass [kg]
% OUTPUT:
% xyzdot - [9,1] Current State Derivatives 

% Load parameters
muE = params.muE;
rE = params.rE;
wE = params.wE;
Cd = params.Cd;
A = params.A;
m = params.m;

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
rho0 = xyz(7); 
h0 = xyz(8);
H = xyz(9);

% Densities
h = r-rE;
rho = rho0.*exp((h0-h)/H);
%drhodh = -rho/H;

% State Derivatives
xyzdot(1:3) = vs;
xyzdot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;
%xyzdot(7:9) = zeros(3,1);

end