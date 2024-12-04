function Xn = stateNoise(X,sigmaState,sigmaMod,nSats,nSatSt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
modStart = nSats*nSatSt + 1;

Xn = zeros(size(X));

for i = 1:nSats
    ind = (i-1)*nSatSt+1:i*nSatSt;
    Xn(:,ind) = ...
        X(:,ind) + (sigmaState .* randn(N,nSatSt));
end

Xn(:,modStart:end) = X(:,modStart:end) + (sigmaMod.*randn(N,length(sigmaMod)));

end