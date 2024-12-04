function [a_drag,varargout] = accelF_drag_exp_6s1_2m(state,params)
% [r_x r_y r_z v_x v_y v_z rho_0 H]
CD = params.C_D;
A = params.A;
m = params.m;
h_0 = params.h_0;
wE = params.wE;
rE = params.rE;

% Angular Velocity
wA = [0;0;wE];

% Extract state
rs = state(1:3); r = sqrt(rs.'*rs);
vs = state(4:6)';
rho_0 = state(7); 
H = state(8);

% Relative velocity
vA = cross(wA,rs);
vRs = vs - vA;
vR = sqrt(vRs.'*vRs);

% Calculate density
h = r - rE;
rho = rho_0*exp((h_0-h)/H);

% Drag acceleration
Dconst = -(CD*A*vR)/(2*m);
a_drag = Dconst*rho*vRs;

if nargout == 4
    rhat = rs./r;
    drhodh = -rho/H;
    vRhat = vRs./vR;
    % Drag acceleration partial derivatives
    dadr = Dconst*drhodh*(vRs.*rhat.');
    dadv = Dconst*rho*(I + vRhat.*vRhat.');
    dadS = [];
    dadM = Dconst*vRs*[exp((h_0-h)/H), -drhodh*((h-h_0)/H)];
    varargout{2} = dadM;
    varargout{1} = [dadr;dadv;dadS];
end


end