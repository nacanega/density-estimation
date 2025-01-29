function xyzPhidot = STM(t, xyzPhi, params)

% Load parameters
muE = params.muE;
rE = params.rE;
wE = params.wE;
Cd = params.Cd;
A = params.A;
m = params.m;
n = params.n;
H_m = params.H_m;
H_M = params.H_M;
h_low = params.h_low;
JD_epoch = params.JDe;
JD = JD_epoch + t / 86400;


% Angular velocity vector
wA = [0; 0; wE];

% Preallocate
xyzPhidot = zeros(size(xyzPhi));

% Extract states
rs = xyzPhi(1:3); % Position
vs = xyzPhi(4:6); % Velocity
rho_m_low = xyzPhi(7); 
rho_M_low = xyzPhi(8); 
Phi = reshape(xyzPhi(9:end), [8, 8]); % STM

% Compute derived quantities
r = norm(rs); 
vR = norm(vs - cross(wA, rs)); 
vRhat = (vs - cross(wA, rs)) / vR; 

h = r - rE;

% Compute density using the exponential model
Y = pqw_to_geo(xyzPhi(1:6), params);
[chi_term,Ub] = Ubfun(Y(1:3),JD,params);
rho_min = rho_m_low .* exp((h_low - h) ./ H_m);
rho_max = rho_M_low .* exp((h_low - h) ./ H_M);
rho = rho_min + (rho_max-rho_min)*chi_term;

% State Derivatives
xyzPhidot(1:3) = vs; 
xyzPhidot(4:6) = -(muE / r^3) * rs - (0.5 * rho * Cd * A / m) * vR^2 * vRhat; 
xyzPhidot(7) = 0; 
xyzPhidot(8) = 0; 

% Gravity contribution
I3 = eye(3);
rhat = rs / r;
dadrG = -(muE / r^3) * (I3 - 3 * (rhat * rhat.')); 

% Drag contribution
x = rs(1); y = rs(2); z = rs(3); vx = vs(1); vy = vs(2); vz = vs(3); Ubx = Ub(1); Uby = Ub(2); Ubz = Ub(3);
dadrD = [(A^2*Cd^2*vR^2*vx*((rho_min*x)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 - (x*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vx*((rho_min*y)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 - (y*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vx*((rho_min*z)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 - (z*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2);
         (A^2*Cd^2*vR^2*vy*((rho_min*x)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 - (x*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vy*((rho_min*y)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 - (y*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vy*((rho_min*z)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 - (z*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2); 
         (A^2*Cd^2*vR^2*vz*((rho_min*x)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 - (x*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vz*((rho_min*y)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 - (y*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2), (A^2*Cd^2*vR^2*vz*((rho_min*z)/(H_m*r) + ((1/2^(n/2 + 1))*n*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2 - 1)*(rho_min - rho_max)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 - (z*((Ubx*x + Uby*y + Ubz*z)/r + 1)^(n/2)*(H_M*rho_min - H_m*rho_max))/(2^(n/2)*H_M*H_m*r)))/(4*m^2)];
dadvD = [-(A*Cd*rho*(2*vx^2 + vy^2 + vz^2))/(2*m*vR), -(A*Cd*rho*vx*vy)/(2*m*vR), -(A*Cd*rho*vx*vz)/(2*m*vR);
         -(A*Cd*rho*vx*vy)/(2*m*vR), -(A*Cd*rho*(vx^2 + 2*vy^2 + vz^2))/(2*m*vR), -(A*Cd*rho*vy*vz)/(2*m*vR);
        -(A*Cd*rho*vx*vz)/(2*m*vR), -(A*Cd*rho*vy*vz)/(2*m*vR), -(A*Cd*rho*(vx^2 + vy^2 + 2*vz^2))/(2*m*vR)];
dadMD = [-(A*Cd*vR*vx*(exp((h_low - r + rE)/H_m) - exp(((h_low - r + rE)/H_m))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2)))/(2*m), -(A*Cd*vR*vx*exp(((h_low - r + rE)/H_M))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2))/(2*m); 
         -(A*Cd*vR*vy*(exp((h_low - r + rE)/H_m) - exp(((h_low - r + rE)/H_m))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2)))/(2*m), -(A*Cd*vR*vy*exp(((h_low - r + rE)/H_M))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2))/(2*m);
         -(A*Cd*vR*vz*(exp((h_low - r + rE)/H_m) - exp(((h_low - r + rE)/H_m))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2)))/(2*m), -(A*Cd*vR*vz*exp(((h_low - r + rE)/H_M))*((Ubx*x + Uby*y + Ubz*z)/(2*r) + sym(1/2))^(n/2))/(2*m)];

% Assemble F matrix
F = [zeros(3), eye(3), zeros(3, 1), zeros(3, 1);
     dadrG + dadrD, dadvD, dadMD;   
     zeros(1, 3), zeros(1, 3), 0, 0;             
     zeros(1, 3), zeros(1, 3), 0, 0];              

% STM derivative
Phidot = F * Phi;

% Flatten and store
xyzPhidot(9:end) = reshape(Phidot, [], 1);

end
