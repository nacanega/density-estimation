function F = F_nodrag(xyz)
%F_nodrag calculates the F (or A) matrix for the simple orbit model with no drag
% and uniform density spherical Earth
% INPUTS:
%   xyz - [6,1(or n)] State vector of positions and velocities [km and km/s]
% OUTPUT:
%     F - [6,6,1(or n))] F matrix for a given state vector

muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]

% Set the trivial portions of F
dvdr = zeros(3);
dvdv = eye(3);
dadv = zeros(3);

if isvector(xyz)  
    % Simple, single State Update Matrix output
    x = xyz(1); y = xyz(2); z = xyz(3);
elseif ismatrix(xyz)
    % More complicated, multiple State Update Matrix output
    [~,N] = size(xyz);
    x = reshape(xyz(1,:),1,1,[]);
    y = reshape(xyz(2,:),1,1,[]);
    z = reshape(xyz(3,:),1,1,[]);
    pages = ones(1,1,N);
    dvdr = dvdr.*pages;
    dvdv = dvdv.*pages;
    dadv = dadv.*pages;
else
    % This is an error
    error("Size:notVectorMatrix","Input must be vector or matrix")
end

r = vecnorm([x;y;z]);
dadr = [ ...
    -muE./r.^3 + (3*muE*x.*x)./r.^5, (3*muE*x.*y)./r.^5, (3*muE*x.*z)./r.^5; ...
    (3*muE*x.*y)./r.^5, -muE./r.^3 + (3*muE*y.*y)./r.^5, (3*muE*y.*z)./r.^5; ...
    (3*muE*x.*z)./r.^5, (3*muE*y.*z)./r.^5, -muE./r.^3 + (3*muE*z.*z)./r.^5 ...
];

% Output constructed matrix
F = [dvdr dvdv; dadr dadv];

end