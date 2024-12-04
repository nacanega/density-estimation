function statedot = stateEQ_drag_Ns1_Nm(~,state,params)
%
%
%

gravFunc = params.gravFunc;
dragFunc = params.dragFunc;

statedot = zeros(size(state));
statedot(1:3) = state(4:6);

a_grav = gravFunc(state,params);
a_drag = dragFunc(state,params);
statedot(4:6) = a_grav + a_drag;

end