function [rho, dadr, dadv, dadM] = hpDensity(t,z,params)
%hpDensity Returns the density and partial derivatives of drag acceleration
% for the Harris-Priester density atmospheric model
%   INPUTS:
%       t - Current time [s]
%       z - (8,1) State vector [r_x r_y r_z v_x v_y v_z rho_min rho_max]'
%  params - Structure of constant parameters and function handles
%           .    muE - Earth gravitational parameter [km^3/s^2]
%           .     rE - Earth radius [km]
%           .     wE - Earth angular velocity [rad/s]
%           .     Cd - Drag Coefficient
%           .      A - Area [km^2]
%           .      m - Mass [kg]
%           . rhoFun - Density Function Handle
%           .      n - Latitudinal Density Variation Exponent
%           .  UbFun - Diurnal Bulge Unit Vector
%   OUTPUT:
%     rho - Harris Priester Density [kg/km^3]
%    dadr - (3,3) Drag Acceleration partial derivatives with respect to
%                 position
%    dadv - (3,3) Drag Acceleration partial derivatives with respect to
%                 velocity
%    dadM - (3,2) Drag Acceleration partial derivatives with respect to
%                 density model parameters (rho_min, rho_max)
%
% From Section 4.5.6 Modified Harris-Priester Atomspheric Model
%
% Long, A.C., Cappellari, Jr. J. O., Velez, C. E., and Fuchs, A. J. (1989) 
%     Goddard Trajectory Determination System (GTDS) Mathematical 
%     Theory (Rev. 1, pp.4-57..4-64) CRC and NASA. 

persistent bands 
persistent heights
persistent min_density
persistent max_density
if isempty("bands")
    bands = [100 120:10:300, 320:20:800, 840:40:1000, Inf];
    heights = bands(1:end-1);
    min_density = [497400, 24900, 8377, 3899, 2122, 1263, 800.8, 528.3, 361.7, 255.7, 183.9, 134.1, 99.49, 74.88, ...
        57.09, 44.03, 34.3, 26.97, 21.39, 17.08, 10.99, 7.214, 4.824, 3.274, 2.249, 1.558, 1.091, ...
        0.7701, 0.5474, 0.3915, 0.2813, 0.2042, 0.1488, 0.1092, 0.0807, 0.06012, 0.04519, 0.0343, ...
        0.02632, 0.02043, 0.01607, 0.01281, 0.01036, 0.008496, 0.007069, 0.00468, 0.0032, 0.00221, ...
        0.00156, 0.00115]/1e3;
    max_density = [497400, 24900, 8710, 4059, 2215, 1344, 875.8, 601, 429.7, 316.2, 239.6, 185.3, 145.5, 115.7, ...
        93.08, 75.55, 61.82, 50.95, 42.26, 35.26, 25.11, 18.19, 13.37, 9.955, 7.492, 5.384, 4.355, ...
        3.362, 2.612, 2.042, 1.605, 1.267, 1.005, 0.7997, 0.639, 0.5123, 0.4121, 0.325, 0.2591, ...
        0.2185, 0.1779, 0.1452, 0.119, 0.09776, 0.08059, 0.05741, 0.0421, 0.0313, 0.0236, 0.0181]/1e3;
end

rE = params.rE;
Cd = params.Cd;
A = params.A;
m = params.m;
n = params.n;
UbFun = params.UbFun; % Input: (t,[rs;vs],params); Output: [Ubx;Uby;Ubz]

% Position
rs = z(1:3);
x = rs(1); y = rs(2); z = rs(3);
r = sqrt(rs.'*rs);
rhat = rs./r;

% Velocity
vs = z(4:6);
v = sqrt(vs.'*vs);
vx = vs(1); vy = vs(2); vz = vs(3);
vhat = vs./v;

% Density Parameters
rho_min = z(7);
rho_max = z(8);

% Altitude
h = r - rE;

% Reference altitudes
band = discretize(h,bands);
h_low = heights(band);       % Lower altitude bound
h_high = heights(band+1);    % Upper altitude bound

% Reference densities
rho_m_low = min_density(band);       % Min density at h_i
rho_m_high = min_density(band+1);    % Min density at h_{i+1}
rho_M_low = max_density(band);       % Max density at h_i
rho_M_high = max_density(band+1);    % Max density at h_{i+1}

% Compute scale heights H_m and H_M
H_m = (h_low - h_high) ./ log(rho_m_high ./ rho_m_low);
H_M = (h_low - h_high) ./ log(rho_M_high ./ rho_M_low);

% Compute density
rho = [];

% Compute Ub
Ub = UbFun(t,z,params);
Ubx = Ub(1);
Uby = Ub(2); 
Ubz = Ub(3);
% Compute partial derivatives of drag acceleration
drconst = (A*Cd*v)/(2*m); Ubconst = (Ubx*x + Uby*y + Ubz*z); 
dadr = drconst*[(vx*((rho_m_low*x*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 + (x*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r))),vx*((rho_m_low*y*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 + (y*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r)), vx*((rho_m_low*z*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 + (z*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r));
                 vy*((rho_m_low*x*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 + (x*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r)), vy*((rho_m_low*y*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 + (y*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r)), vy*((rho_m_low*z*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 + (z*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r));
                 vz*((rho_m_low*x*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubx*y^2 - Uby*x*y + Ubx*z^2 - Ubz*x*z))/r^3 + (x*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r)), vz*((rho_m_low*y*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Uby*x^2 - Ubx*y*x + Uby*z^2 - Ubz*y*z))/r^3 + (y*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r)), vz*((rho_m_low*z*exp(((h_low + rE - r)/H_m)))/(H_m*r) - ((1/2^(n/2 + 1))*n*(rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2 - 1)*(Ubz*x^2 - Ubx*z*x + Ubz*y^2 - Uby*z*y))/r^3 + (z*(H_m*rho_M_low*exp(((h_low + rE - r)/H_M)) - H_M*rho_m_low*exp(((h_low + rE - r)/H_m)))*(Ubconst/r + 1)^(n/2))/(2^(n/2)*H_M*H_m*r))];
dadv = (drconst/v)*[- (((rho_M_low*exp(((h_low + rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m)))*v) - (vx^2*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), -(vx*vy*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), -(vx*vz*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v);
                    -(vx*vy*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), - (((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m)))*v) - (vy^2*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), -(vy*vz*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v);
                    -(vx*vz*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), -(vy*vz*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v), - (((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m)))*v) - (vz^2*((rho_M_low*exp(((h_low+ rE - r)/H_M)) - rho_m_low*exp(((h_low+ rE - r)/H_m)))*((Ubconst)/(2*r) + sym(1/2))^(n/2) + rho_m_low*exp(((h_low+ rE - r)/H_m))))/(v)];
dadM = [];

end