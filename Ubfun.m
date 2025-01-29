function [chi_term,UB] = Ubfun(r,JD,params)
n = params.n;
tilt = params.tilt; % Earth's axial tilt in degrees 
lambda_deg = params.lambda; % Lag angle in degrees 
lambda_rad = deg2rad(lambda_deg);
JD_epoch = params.JDe; % Julian Date at epoch

% Fractional year (gamma) in radians
gamma = 2 * pi * (JD - JD_epoch) / 365.25;

% Declination of the Sun (delta_s) in radians
delta_s = asin(sin(deg2rad(tilt)) * sin(gamma));

% Right ascension of the Sun (alpha_s) in radians
alpha_s = atan2(cos(deg2rad(tilt)) * sin(gamma), cos(gamma));

% Apex right ascension (alpha_B), normalized to [0, 2*pi]
alpha_B = mod(alpha_s + lambda_rad, 2*pi);

% Compute U_B (unit vector toward the apex of the diurnal bulge)
UBx = cos(delta_s) * cos(alpha_B);
UBy = cos(delta_s) * sin(alpha_B);
UBz = sin(delta_s);

UB = [UBx; UBy; UBz]; % Unit vector in geocentric equatorial coordinates

chi_term = (.5 + (dot(r,UB)./(2*norm(r))))^(n/2);

end
