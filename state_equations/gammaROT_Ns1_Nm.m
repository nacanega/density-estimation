function [Gamma,varagout] = gammaROT_Ns1_Nm(dt,state)
%gammaROT_7s returns the Gamma and Rotation matrices for a 7 state system
nargoutchk(1,2)

N = length(state);
I = eye(3);
Gamma = zeros(N,3);
Gamma(1:6,:) = dt*[(dt/2)*I;I];

if nargout == 2
    rs = state(1:3); r = sqrt(rs.'*rs); rhat = rs./r;
    vs = state(4:6)'; v = sqrt(vs.'*vs); vhat = vs./v;
    nhat = cross(-rhat,vhat);
    varagout{1} = [-rhat,vhat,nhat];
end

end