function statedot = stateSTM_drag_Ns1_Nm(~,state,params)
%
%
%

gravFunc = params.gravFunc;
dragFunc = params.dragFunc;
nSatSt = params.nSatStates;
nModSt = params.nModStates;

statedot = zeros(size(state));
statedot(1:3) = state(4:6);

[a_grav,dadrvSG,dadMG] = gravFunc(state,params);
[a_drag,dadrvSD,dadMD] = dragFunc(state,params);

statedot(4:6) = a_grav + a_drag;
 
dvdrvS = zeros(3,nSatSt);
dvdrvS(:,4:6) = eye(3);
dvdM = zeros(3,nModSt);

dadrvS = dadrvSG + dadrvSD;
dadM = dadMG + dadMD;

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