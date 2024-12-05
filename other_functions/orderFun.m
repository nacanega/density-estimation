function orders = orderFun(X0,nSats,nSatSt,nModSt)
%ORDERFUN Summary of this function goes here
%   Detailed explanation goes here
N = nSats*nSatSt + nModSt;

orders = zeros(N,1);

for i = 1:nSats
    j = (i-1)*nSatSt;
    orders((1:3)+j) = 10.^max(floor(log10(norm(X0((1:3)+j)))));
    orders((4:6)+j) = 10.^max(floor(log10(norm(X0((4:6)+j)))));
    orders((7:nSatSt)+j)  = 10.^floor(log10(abs(X0((7:nSatSt)+j))));
end

orders(end-(nModSt-1):end) = ...
    10.^(fix(log10(abs(X0(end-(nModSt-1):end)))-1)+1);

end

