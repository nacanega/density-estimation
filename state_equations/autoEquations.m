function [aGFun,aDFun,gFun,seqFun,stmFun] = autoEquations(satStates,modStates,models,nStates)
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

dsuffix = sprintf("_%ds1_%dm",m,n);
gsuffix = sprintf("_%ds1_Nm",m);

if (nStates - n) >= 2*m
    % Block functions
    gamMod = "_NsN_Nm";
    evalGam = true;
else % Individual functions
    gamMod = "_Ns1_Nm";
    evalGam = false;
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
    case {"exp","Exp","exponential","Exponential"}
        rhoModelName = "_exp";
        dragType =  "_drag";
    case {"hp","HP","harris-priester","Harris-Priester"}
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
gAccString = gAccString + gModelName + gsuffix;
dAccString = dAccString + rhoModelName + dsuffix;

aGFun = str2func(gAccString);
aDFun = str2func(dAccString);
gFun = str2func(gammaString);
if evalGam
   gFun = @(dt,X) gFun(dt,X,m,n);
end
seqFun = str2func(seqString);
stmFun = str2func(stmString);

end