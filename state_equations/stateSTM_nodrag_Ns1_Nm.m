function statedot = stateSTM_nodrag_Ns1_Nm(~,state,params)
%
%
%

gravFunc = params.gravFunc;
nSatSt = params.nSatStates;
nModSt = params.nModStates;

statedot = zeros(size(state));
statedot(1:3) = state(4:6);

[a_grav,dadrvS,dadM] = gravFunc(state,params);

statedot(4:6) = a_grav;
 
dvdrvS = zeros(3,nSatSt);
dvdrvS(:,4:6) = eye(3);
dvdM = zeros(3,nModSt);

ZS = zeros(nSatSt-6,nSatSt);
ZM = zeros(nModSt,nSatSt);

F = sparse([ ...
    dvdrvS,dvdM; ...
    dadrvS,dadM;
    ZS; ZM ]);

N = length(F);
Phi = reshape(state(N+1:end),N,N);
statedot(N+1:end) = reshape(F*Phi,[],1);

end