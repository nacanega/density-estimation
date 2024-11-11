function xyzPhidot = stateSTM_phasedrag_9s8p(~,xyzPhi,params)
%stateSTM_phasedrag_9s8p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyzPhi - [9+9*9,1] Current State and STM
%    params - Structure of parameters
%           .    muE - Earth gravitational parameter [km^3/s^2]
%           .     rE - Earth radius [km]
%           .     wE - Earth angular velocity [rad/s]
%           . rhoFun - Density function handle
%           .     Cd - Drag Coefficient
%           .      A - Area [km^2]
%           .      m - Mass [kg]
%           . lambda - Sun parameter [?]
% OUTPUT:
% xyzPhidot - [9+9*9,1] Current State and STM Derivatives 

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
xyzPhidot = zeros(size(xyzPhi));

% Positions
rx = xyzPhi(1); ry = xyzPhi(2); rz = xyzPhi(3);
rs = [rx;ry;rz];
r = norm(rs);
rhat = rs./r;

% Velocities
vs = xyzPhi(4:6);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = norm(vRs);
vRhat = vRs./vR;

% Density
[rho,drhodh] = rhoFun(~,r,lambda);

% State Derivatives
xyzPhidot(1:3) = vs;
xyzPhidot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;

% Phi = reshape(xyzPhi(7:end),[6 6])
Z = zeros(3);
I = eye(3);

% Velocity
dvdr = Z;
dvdv = I;

% Gravity contribution
dadvG = Z;
dadrG = -(muE/r^3)*(I - 3*rhat.*rhat.');

% Drag contribution
Dconst = -((Cd*A*vR)/(2*m));
dadrD = Dconst*drhodh*vRs.*rhat.';
dadvD = Dconst*rho*(I + vRhat.*vRhat.');

% Sum acceleration components
dadr = dadrG + dadrD;
dadv = dadvG + dadvD;

% Assemble F matrix
F = [dvdr,dvdv;...
     dadr,dadv];

% State derivatives
xyzPhidot(7:end) = reshape(F,size(xyzPhi(7:end)));

end