function modPs = densityParams(r,varargin)
%expDensity Implements the basic exponential atmospheric model from Vallado
% which incorporates several bands
%   INPUTS:
%         r - [m,1] radius [km]
%  varargin - placeholder to allow for interchangeability with more complex
%             atmospheric density models
%   OUTPUT:
%     modPs - [m,3] vector of model parameters
%      rho0 - [m,1] nominal density [kg/km^3]
%        h0 - [m,1] base altitude [km]
%         H - [m,1] scale height [km]
%
% From Table 8-4. Exponential Atmospheric Model 
%
% Vallado, D. A., and Wertz, J.R. (2022) Model Atmospheres. 
%   In Fundamentals of Astrodynamics and Applications 
%   (5th ed., pp. 568–578). Microcosm Press.

rE = 6378.1363; % [km]

% Parameters for each band
bands = [0 25 30:10:150 180 200:50:450 500:100:1000 Inf];
h0s = bands(1:end-1);
rho0s = [1.225e+9 3.899e+7 1.774e+7 3.972e+6 1.057e+6 3.206e+5 ...
         8.770e+4 1.905e+4 3.396e+3 5.297e+2 9.661e+1 2.438e+1 ...
         8.484e+0 3.845e+0 2.070e+0 5.464e-1 2.789e-1 7.248e-2 ...
         2.418e-2 9.518e-3 3.725e-3 1.585e-3 6.967e-4 1.454e-4 ...
         3.614e-5 1.170e-5 5.245e-6 3.019e-6];
Hs = [7.249  6.349  6.682  7.554  8.382  7.714  6.549   5.799  5.382 ...
      5.877  7.263  9.473 12.636 16.149 22.523 29.740  37.105 45.546 ...
     53.628 53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00];

% Determine which band to use
band = discretize(r-rE,bands);

% Lookup values for that band
rho0 = rho0s(band); % Nominal density
h0 = h0s(band);     % Base altitude
H = Hs(band);       % Scale height

modPs = [rho0 h0 H]; % Output Model Parameters

end