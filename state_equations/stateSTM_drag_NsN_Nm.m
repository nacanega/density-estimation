function statedot = stateSTM_drag_NsN_Nm(~,state,params)

gravFunc = params.gravFunc;
dragFunc = params.dragFunc;
nSatSt = params.nSatStates;
nModSt = params.nModStates;
nSats = params.nSats;

N = nSats*nSatSt + nModSt;
statedot = zeros(size(state));
istatedot = zeros(nSatSt,1);

dvaSdrvS = zerosCell([nSats,1],[nSatSt,nSatSt]);

dvdrvS = zeros(3,nSatSt); dvdrvS(:,4:6) = eye(3);
dvaSdM = zeros(nSats*nSatSt,nModSt);
dMdrvS = zeros(nModSt,N);
ZS = zeros(nSatSt-6,nSatSt);

for i = nSats:-1:1
    n1 = nSatSt*(i-1) + 1; n2 = i*nSatSt;
    istate = [state(n1:n2); state(end-(nModSt-1):end)];
    [a_grav,dadrvSG,dadMG] = gravFunc(istate,params);
    [a_drag,dadrvSD,dadMD] = dragFunc(istate,params);
    a_net = a_grav + a_drag;
    dadrvS = dadrvSG + dadrvSD;
    dadM = dadMG + dadMD;
    istatedot(1:6) = [ ...
        state(n1+3:n2-(nSatSt-6));
        a_net
    ];
    statedot(n1:n2) = istatedot;
    dvaSdrvS{i} = [dvdrvS; dadrvS; ZS];
    dvaSdM(n1+3:n1+5,:) = dadM;
end

F = sparse([blkdiag(dvaSdrvS{:}),dvaSdM;dMdrvS]);

Phi = reshape(state(N+1:end),N,N);
statedot(N+1:end) = reshape(F*Phi,[],1);

end