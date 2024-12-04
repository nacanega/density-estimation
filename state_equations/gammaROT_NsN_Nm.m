function [Gamma,varargout] = gammaROT_NsN_Nm(dt,state,nSatStates,nModStates)
%gammaROT_7s returns the Gamma and Rotation matrices for a 7 state system
nargoutchk(1,2)

N = length(state);
nSats = (N - nModStates) / nSatStates;

Gamma = zeros(N,3*nSats);

if nargout == 2
    gamma = zerosCell([nSats,1],[nSatStates,3]);
    rotms = zerosCell([nSats,1],[3,3]);
    for i = nSats:-1:1
        iState = state(nSatStates*(i-1)+1:nSatStates*i);
        [gamma{i},rotms{i}] = gammaROT_Ns1_Nm(dt,iState);
    end
    Gamma(1:end-mModStates,:) = sparse(blkdiag(gamma{:}));
    varargout = sparse(blkdiag(rotms{:}));
else
    gamma = zerosCell([nSats,1],[nSatStates,3]);
    for i = nSats:-1:1
        iState = state(nSatStates*(i-1)+1:nSatStates*i);
        gamma{i} = gammaROT_Ns1_Nm(dt,iState);
    end
    Gamma(1:end-nModStates,:) = sparse(blkdiag(gamma{:}));
end

end