function xyzdot = stateEQ_nodrag(~,xyz,params)
%stateEQ_NLPhi is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyz - [6,1] Current State
% params - Structure of parameters
%        .    muE - Earth gravitational parameter [km^3/s^2]
% OUTPUT:
% xyzdot - [6,1] Current State Derivatives 

% Load parameters
muE = params.muE; 

% Preallocate
xyzdot = zeros(size(xyz));

% Positions
rs = xyz(1:3); vs = xyz(4:6);
r = sqrt(rs.'*rs);

% State Derivatives
xyzdot(1:3) = vs;
xyzdot(4:6) = -(muE/r^3)*rs;

end