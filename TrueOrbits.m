function Xdot = TrueOrbits(t,X,params)


% Load parameter
muE = params.muE;
rE = params.rE;
wE = params.wE;
Cd = params.Cd;
A = params.A;
m = params.m;
JD_epoch = params.JDe;
JD = JD_epoch + t / 86400;

% Angular velocity vector
wA = [0;0;wE];

% Preallocate
Xdot = zeros(size(X));

% Positons
rs = X(1:3);
r = sqrt(rs.'*rs);

% Velocities
vs = X(4:6);
vAs = cross(wA,rs);
vRs = vs-vAs;
vR = norm(vRs);

% Densities
h = r-rE;
Y = pqw_to_geo(X(1:6), params);
[chi_term,~] = Ubfun(Y(1:3),JD,params);
[rho,~,~,~,~] = hpDensity(h,chi_term);

% State Derivatives
Xdot(1:3) = vs;
Xdot(4:6) = -(muE/r^3)*rs - (0.5*rho*Cd*A/m)*vR*vRs; % 
end