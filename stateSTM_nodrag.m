function xyzPhidot = stateSTM_nodrag(~,xyzPhi)
%stateEQ_NLPhi is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyzPhi - [6+6*6,1] Current State and STM
% OUTPUT:
% xyzPhidot - [6+6*6,1] Current State and STM Derivatives 

muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]

% Temporarily define parameters here
% FIXME Add parameters input

xyzPhidot = zeros(size(xyzPhi));
rx = xyzPhi(1); ry = xyzPhi(2); rz = xyzPhi(3);
rs = [rx;ry;rz];
r = norm(rs);

% State Derivatives
xyzPhidot(1:3) = xyzPhi(4:6);
xyzPhidot(4:6) = -(muE/r^3)*rs;

% Phi = reshape(xyzPhi(7:end),[6 6])
dvdr = zeros(3);
dvdv = eye(3);
dadv = zeros(3);
dadr = [ ...
    -muE/r^3 + (3*muE*rx*rx)/r^5, (3*muE*rx*ry)/r^5, (3*muE*rx*rz)/r^5; ...
    (3*muE*rx*ry)/r^5, -muE/r^3 + (3*muE*ry*ry)/r^5, (3*muE*ry*rz)/r^5; ...
    (3*muE*rx*rz)/r^5, (3*muE*ry*rz)/r^5, -muE/r^3 + (3*muE*rz*rz)/r^5 ...
];
F = [dvdr dvdv; dadr dadv];
xyzPhidot(7:end) = reshape(F,size(xyzPhi(7:end)));

end