function [a_grav,varargout] = accelF_grav_6s1_Nm(state,params)
% [r_x r_y r_z v_x v_y v_z C_D rho_0 h_0 H]
muE = params.muE;

% Extract state
rs = state(1:3); r = sqrt(rs.'*rs);

% Gravity acceleration
Gconst = -muE/r^3;
a_grav = Gconst*rs;

if nargout == 4
    rhat = rs./r;
    n = length(state) - 6;
    % Gravity acceleration partial derivatives
    dadr = Gconst*(I - 3*rhat.*rhat.');
    dadv = zeros(3);
    dadS = [];
    dadM = zeros(3,n);
    varargout{2} = dadM;
    varargout{1} = [dadr;dadv;dadS];
end

end