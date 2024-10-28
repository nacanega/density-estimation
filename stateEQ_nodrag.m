function xyzdot = stateEQ_nodrag(~,xyz)
%stateEQ_NLPhi is a function defining the state derivatives as a function
% of the existing state
% INPUTS:
%    xyz - [6,1] Current State
% OUTPUT:
% xyzdot - [6,1] Current State Derivatives 

muE = 398600.4415; % Earth gravitational parameter [km^3/s^2]
% Temporarily define parameters here
% FIXME Add parameters input

xyzdot = zeros(size(xyz));
rs = xyz(1:3); vs = xyz(4:6);
r = sqrt(rs.'*rs);

% State Derivatives
xyzdot(1:3) = vs;
xyzdot(4:6) = -(muE/r^3)*rs;

end