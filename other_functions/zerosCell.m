function output = zerosCell(cellSize,matSize)
%zerosCell Preallocates a cell array (k-rows, 1 column) with
% INPUTS:
% cellSize - [m1, m2, ... mk] Vector of cell array dimension lengths
%  matSize - [n1, n2, ... nk] Vector of zero array dimension lengths
% OUTPUT:
% zerosCell - [m1, m2, ... mk] Cell array of [n1, n2, ... nk] zero matrices

output = cell(cellSize);
output(:) = {zeros(matSize)};

end