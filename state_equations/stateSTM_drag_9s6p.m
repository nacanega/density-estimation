function xyzPhidot = stateSTM_drag_9s6p(~,xyzPhi,params)
%stateSTM_drag_9s6p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyzPhi - [9+9*9,1] Current State and STM
%    params - Structure of parameters
%           .    muE - Earth gravitational parameter [km^3/s^2]
%           .     rE - Earth radius [km]
%           .     wE - Earth angular velocity [rad/s]
%           .     Cd - Drag Coefficient
%           .      A - Area [km^2]
%           .      m - Mass [kg]
% OUTPUT:
% xyzPhidot - [9+9*9,1] Current State and STM Derivatives 

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
xyzPhidot = zeros(size(xyzPhi));

% Density parameters
rho0 = xyz(1); 
h0 = xyz(2);
H = xyz(3);

% Positions
rx = xyzPhi(4); ry = xyzPhi(5); rz = xyzPhi(6);
rs = [rx;ry;rz];
r = norm(rs);
rhat = rs./r;

% Velocities
vs = xyzPhi(7:9);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = norm(vRs);
vRhat = vRs./vR;

% Density
h = r-rE;
rho = rho0.*exp((h0-h)/H);
drhodh = -rho/H;

% State Derivatives
xyzPhidot(1:3) = zeros(3,1);
xyzPhidot(4:6) = vs;
xyzPhidot(7:9) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;

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

% Model Contribution
dadM = Dconst*rho*vRs./[rho0,H,H*(H/(h-h0))];

% Sum acceleration components
dadr = dadrG + dadrD;
dadv = dadvG + dadvD;

% Assemble F matrix
F = [   Z,   Z,   Z;...
        Z,dvdr,dvdv;...
     dadM,dadr,dadv];

% State derivatives
xyzPhidot(10:end) = reshape(F,size(xyzPhi(10:end)));

end