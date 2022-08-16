%--- Description ---%
%
% Filename: chebmat.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the 1D Chebyshev matrix
%
% Inputs:
% grid - a column vector of points
% k - the desired number of polynomials to use
%
% Output:
% A - the matrix of the first k Chebyshev polynomials evaluated on the grid

function A = chebmat(grid,k)
A = zeros(length(grid),k);

for i = 1:k
    if i == 1
        A(:,i) = 1;
    else
        A(:,i) = sqrt(2).*cos((i-1).*acos(grid));
    end
end
end