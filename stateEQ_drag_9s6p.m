function xyzdot = stateEQ_drag_6s7p(~,xyz,params)
%stateEQ_drag_6s7p is a function defining the state derivatives as a function
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

% Density parameters
rho0 = xyz(1); 
h0 = xyz(2);
H = xyz(3);

% Positons
rs = xyz(4:6);
r = sqrt(rs.'*rs);

% Velocities
vs = xyz(7:9);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = norm(vRs);

% Densities
h = r-rE;
rho = rho0.*exp((h0-h)/H);
%drhodh = -rho/H;

% State Derivatives
xyzdot(1:3) = zeros(3,1);
xyzdot(4:6) = vs;
xyzdot(7:9) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;

end