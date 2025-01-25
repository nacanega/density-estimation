function H = HmatPos(nSat,nSatState,nModState)
%HmatPos HmatPos returns the observation (H) matrix for Kalman filtering 
% assuming that the first three values of each satellite state are the 
% positions and that all positions are observed.
%   INPUTS:
%      nSat - Number of satellites
% nSatState - Number of states for individual satellites
% nModState - Number of shared model states
%   OUTPUT:
%         H - Sparse observation matrix

H = zeros(nSat*3,nSat*nSatState+nModState);

for i = 1:nSat
    H(3*(i-1)+(1:3),nSatState*(i-1)+(1:3)) = eye(3);
end

H = sparse(H);
end