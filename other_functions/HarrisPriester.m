function [rho_min, rho_max] = HarrisPriester(h)
    % HarrisPriester: Interpolate atmospheric density based on Harris-Priester model
    % Input:
    % h - Array of altitudes (km)
    % Output:
    % rho_min - Minimum interpolated density (kg/km^3)
    % rho_max - Maximum interpolated density (kg/km^3)

    % Load the table (replace with actual filename if using an Excel file)
    % Table format: [Height (km), Min Density (gm/km^3), Max Density (gm/km^3)]
 % Define the data as three columns
heights = [100, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, ...
        300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, ...
        660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000];
min_density = [497400, 24900, 8377, 3899, 2122, 1263, 800.8, 528.3, 361.7, 255.7, 183.9, 134.1, 99.49, 74.88, ...
        57.09, 44.03, 34.3, 26.97, 21.39, 17.08, 10.99, 7.214, 4.824, 3.274, 2.249, 1.558, 1.091, ...
        0.7701, 0.5474, 0.3915, 0.2813, 0.2042, 0.1488, 0.1092, 0.0807, 0.06012, 0.04519, 0.0343, ...
        0.02632, 0.02043, 0.01607, 0.01281, 0.01036, 0.008496, 0.007069, 0.00468, 0.0032, 0.00221, ...
        0.00156, 0.00115];
max_density = [497400, 24900, 8710, 4059, 2215, 1344, 875.8, 601, 429.7, 316.2, 239.6, 185.3, 145.5, 115.7, ...
        93.08, 75.55, 61.82, 50.95, 42.26, 35.26, 25.11, 18.19, 13.37, 9.955, 7.492, 5.384, 4.355, ...
        3.362, 2.612, 2.042, 1.605, 1.267, 1.005, 0.7997, 0.639, 0.5123, 0.4121, 0.325, 0.2591, ...
        0.2185, 0.1779, 0.1452, 0.119, 0.09776, 0.08059, 0.05741, 0.0421, 0.0313, 0.0236, 0.0181];


 for i = 1:length(h)
        % Find the altitude interval [h_i, h_{i+1}]
        idx = find(heights <= h(i), 1, 'last');
        
        if isempty(idx) || idx == length(heights)
            error('Altitude out of range. Must be between %d and %d km.', min(heights), max(heights));
        end
        
        h_i = heights(idx);       % Lower altitude bound
        h_next = heights(idx+1);  % Upper altitude bound
       
        rho_m_i = min_density(idx);       % Min density at h_i
        rho_m_next = min_density(idx+1); % Min density at h_{i+1}
        rho_M_i = max_density(idx);       % Max density at h_i
        rho_M_next = max_density(idx+1); % Max density at h_{i+1}

        % Compute scale heights H_m and H_M
        H_m = (h_i - h_next) / log(rho_m_next / rho_m_i);
        H_M = (h_i - h_next) / log(rho_M_next / rho_M_i);

        % Exponential interpolation for current altitude
        rho_min(i) = rho_m_i * exp((h_i - h(i)) / H_m);
        rho_max(i) = rho_M_i * exp((h_i - h(i)) / H_M);
 
end
end