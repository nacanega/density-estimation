function zs = measNoise(X,sigmaMeas,nSats,H)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[N,~] = size(X);
[obs,iStates] = size(H);

indexes = H*(1:iStates)';
zs = zeros(N,obs*nSats);

for i = 1:nSats
    ind = indexes + (i-1)*iStates;
    zs(:,(i-1)*obs+1:i*obs) = X(:,ind) + (sigmaMeas' .* randn(N,obs));
end

end