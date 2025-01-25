function zs = measNoise(X,sigmaMeas,nSats,H)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[N,~] = size(X);
[obs,iStates] = size(H);

sigmaMeas = repmat(sigmaMeas,nSats,1);

ind = H*(1:iStates)';
zs = X(:,ind) + (sigmaMeas' .* randn(N,obs));

end