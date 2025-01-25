function [Gamma,varargout] = gammaROT_Ns1_Nm(dt,state)
%gammaROT_7s returns the Gamma and Rotation matrices for a 7 state system
nargoutchk(1,2)

N = length(state);
I = eye(3);
Gamma = zeros(N,3);
Gamma(1:6,:) = dt*[(dt/2)*I;I];

if nargout == 2
    rs = state(1:3); rhat = rs./sqrt(rs.'*rs);
    vs = state(4:6); % vhat = vs./sqrt(vs.'*vs);
    hs = cross(rs,vs); hhat = hs./sqrt(hs.'*hs);
    ns = cross(-hs,-rs); nhat = ns./sqrt(ns.'*ns);
    varargout{1} = [nhat,-hhat,-rhat]; %LVLH
end

end