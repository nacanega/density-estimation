function statedot = stateEQ_nodrag_NsN_Nm(~,state,params)

gravFunc = params.gravFunc;
nSatSt = params.nSatStates;
nModSt = params.nModStates;
nSats = params.nSats;

N = nSats*nSatSt + nModSt;
statedot = zeros(N,1);
istatedot = zeros(nSatSt,1);

for i = nSats:-1:1
    a_grav = gravFunc(state,params);
    istatedot(1:6) = [ ...
        state(nSatSt*(i-1)+4:i*nSatSt);
        a_grav
    ];
    statedot(nSatSt*(i-1)+1:i*nSatSt) = istatedot;
end

end