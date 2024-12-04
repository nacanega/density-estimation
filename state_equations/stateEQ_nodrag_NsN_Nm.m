function statedot = stateEQ_nodrag_NsN_Nm(~,state,params)

gravFunc = params.gravFunc;
nSatSt = params.nSatStates;
nModSt = params.nModStates;
nSats = params.nSats;

statedot = zeros(size(state));
istatedot = zeros(nSatSt,1);

for i = nSats:-1:1
    n1 = nSatSt*(i-1) + 1; n2 = i*nSatSt;
    istate = [state(n1:n2); state(end-(nModSt-1):end)];
    a_grav = gravFunc(istate,params);
    istatedot(1:6) = [ ...
        state(n1+3:n2-(nSatSt-6));
        a_grav
    ];
    statedot(nSatSt*(i-1)+1:i*nSatSt) = istatedot;
end

end