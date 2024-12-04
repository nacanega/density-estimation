function statedot = stateEQ_nodrag_Ns1_Nm(~,state,params)
%
%
%

gravFunc = params.gravFunc;

statedot = zeros(size(state));
statedot(1:3) = state(4:6);

a_grav = gravFunc(state,params);
statedot(4:6) = a_grav;

end