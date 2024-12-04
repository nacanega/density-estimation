function supMat = blockMat(n,subMat,varargin)
% Converts a submatrix into a block diagonal matrix repeated n times

if issparse(subMat)
    subMat = full(subMat);
end

if nargin == 3
    interMats = cell(n+1,1);
    interMats(1:n) = {subMat};
    interMats(n+1) = varargin(1);
else
    interMats = cell(n,1);
    interMats(:) = {subMat};
end

supMat = sparse(blkdiag(interMats{:}));

end