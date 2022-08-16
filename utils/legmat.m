%--- Description ---%
%
% Filename: legmat.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the 1D Legendre matrix
%
% Inputs:
% grid - a column vector of points
% k - the desired number of polynomials to use
%
% Output:
% A - the matrix of the first k Legendre polynomials evaluated on the grid

function A = legmat(grid,k)
A = zeros(length(grid),k);
A(:,1) = 1;
A(:,2) = grid*sqrt(3);
for i = 2:k-1
    A(:,i+1)=(grid.* (2*i - 1).*A(:,i)./sqrt(i-1/2)- (i - 1).*A(:,i-1)./sqrt(i-3/2)).*sqrt(i+1/2)/i;
end
end