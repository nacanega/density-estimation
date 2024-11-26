function xyzPhidot = stateSTM_hp_8s9p(~,xyzPhi,params)
%stateSTM_hp_8s6p is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%         t - Current time [s]
%    xyzPhi - [8+8*8,1] Current State and STM
%    params - Structure of parameters
%           .    muE - Earth gravitational parameter [km^3/s^2]
%           .     rE - Earth radius [km]
%           .     wE - Earth angular velocity [rad/s]
%           .     Cd - Drag Coefficient
%           .      A - Area [km^2]
%           .      m - Mass [kg]
%           . rhoFun - Density Function Handle
%           .      n - Latitudinal Density Variation Exponent
%           .  UbFun - Diurnal Bulge Unit Vector
% OUTPUT:
% xyzPhidot - [8+8*8,1] Current State and STM Derivatives 

% Load parameters
muE = params.muE;
%rE = params.rE;
wE = params.wE;
Cd = params.Cd;
A = params.A; 
m = params.m;
rhoFun = params.rhoFun; % Inputs: (t,[rs;vs;rho_min;rho_max],params)
% n = params.n;
% UbFun = params.UbFun; % Input: (t,[rs;vs],params); Output: [Ubx;Uby;Ubz]

% Angular velocity vector
wA = [0;0;wE];

% Preallocate
xyzPhidot = zeros(size(xyzPhi));

% Positions
rs = xyzPhi(1:3);
r = sqrt(rs.'*rs);

% Velocities
vs = xyzPhi(4:6);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = sqrt(vRs.'*vRs);

% Density parameters
rho_min = xyzPhi(7); 
rho_max = xyzPhi(8);

% Density
[rho, dadrD, dadvD, dadMD] = rhoFun(t,[rs;vs;rho_min;rho_max],params);

% State Derivatives
xyzPhidot(1:3) = vs;
xyzPhidot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs;
%xyzPhidot(7:9) = zeros(3,1);

Phi = reshape(xyzPhi(9:end),8,8);
Z = zeros(2);
Zh = zeros(2,3);
Zv = zeros(3,2);
I = eye(3);

% Velocity
dvdr = zeros(3);
dvdv = I;

% Gravity contribution
dadvG = zeros(3);
dadrG = -(muE/r^3)*(I - 3*rhat.*rhat.');

% % Drag contribution
% Dconst = -((Cd*A*vR)/(2*m));
% dadrD = Dconst*drhodh*vRs.*rhat.';
% dadvD = Dconst*rho*(I + vRhat.*vRhat.');
% 
% % Model Contribution
% dadM = Dconst*rho*vRs./[rho0,H*(H/(h-h0))];

% Sum acceleration components
dadr = dadrG + dadrD;
dadv = dadvG + dadvD;
dadM = dadMD;

% Assemble F matrix
F = [dvdr, dvdv,   Zv; ...
     dadr, dadv, dadM; ...
       Zh,   Zh,    Z];

% State derivatives
xyzPhidot(9:end) = reshape(F*Phi,[],1);

end