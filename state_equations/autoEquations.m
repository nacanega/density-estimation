function [aGFun,aDFun,gFun,seqFun,stmFun] = autoEquations(satStates,modStates,models,combined)
%autoEquations
% INPUTS:
% OUTPUT:
% TODO validate

m = length(satStates);
n = length(modStates);

gAccString = "accelF_grav";
dAccString = "accelF_drag";
gammaString = "gammaROT";
seqString = "stateEQ";
stmString = "stateSTM";

if combined 
    % Block functions
    suffix = sprintf("_%ds1_%dm",m,n);
    gamMod = "_NsN_Nm";
else % Individual functions
    suffix = sprintf("_%dsN_%dm",m,n);
    gamMod = "_Ns1_Nm";
end

switch models.gravity
    case {"","sphere","point","J0","J_0","J1","J_1"}
        gModelName = "";
    case {"J2","J_2","J3","J_3","J4","J_4","J5","J_5","J6","J_6"}
        gModelName = "_JN";
    otherwise
        eid = "Model:unrecognizedGravityPerturbation";
        msg = "Unrecognized gravity perturbation model.";
        error(eid,msg);
end

switch models.density
    case ""
        rhoModelName = "_none";
        dragType =  "_nodrag";
    case {"exp","exponential"}
        rhoModelName = "_exp";
        dragType =  "_drag";
    case {"hp","HP","Harris-Priester"}
        rhoModelName = "_hp";
        dragType = "_drag";
    otherwise
        eid = "Model:unrecognizedDensity";
        msg = "Unrecognized density model.";
        error(eid,msg);
end

gammaString = gammaString + gamMod;
seqString = seqString + dragType + gModelName + gamMod;
stmString = stmString + dragType + gModelName + gamMod;
gAccString = gAccString + gModelName + suffix;
dAccString = dAccString + rhoModelName + suffix;

aGFun = str2func(gAccString);
aDFun = str2func(dAccString);
gFun = str2func(gammaString);
seqFun = str2func(seqString);
stmFun = str2func(stmString);

end