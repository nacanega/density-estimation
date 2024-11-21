function [rho_min, rho_max] = hpDensity(h)
    % HarrisPriester: Interpolate atmospheric density based on Harris-Priester model
    % Input:
    % h - Array of altitudes (km)
    % Output:
    % rho_min - Minimum interpolated density (kg/km^3)
    % rho_max - Maximum interpolated density (kg/km^3)

    % Load the table (replace with actual filename if using an Excel file)
    % Table format: [Height (km), Min Density (gm/km^3), Max Density (gm/km^3)]
 % Define the data as three columns
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

band = discretize(h,bands);

h_low = heights(band);       % Lower altitude bound
h_high = heights(band+1);    % Upper altitude bound

rho_m_low = min_density(band);       % Min density at h_i
rho_m_high = min_density(band+1);    % Min density at h_{i+1}
rho_M_low = max_density(band);       % Max density at h_i
rho_M_high = max_density(band+1);    % Max density at h_{i+1}

% Compute scale heights H_m and H_M
H_m = (h_low - h_high) ./ log(rho_m_high ./ rho_m_low);
H_M = (h_low - h_high) ./ log(rho_M_high ./ rho_M_low);

% Exponential interpolation for current altitude
rho_min = rho_m_low .* exp((h_low - h) ./ H_m);
rho_max = rho_M_low .* exp((h_low - h) ./ H_M);

end