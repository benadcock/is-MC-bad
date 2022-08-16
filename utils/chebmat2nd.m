%--- Description ---%
%
% Filename: chebmat2nd.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the 1D Chebyshev (2nd kind) matrix
%
% Inputs:
% grid - a column vector of points
% k - the desired number of polynomials to use
%
% Output:
% A - the matrix of the first k Chebyshev polynomials (of the 2nd kind) evaluated on the grid

function A = chebmat2nd(grid,k)
A = zeros(length(grid),k);

grid (grid == 0) = 1e-15; % if the grid includes zero, add a small number to avoid division by zero

for i = 1:k
    if i == 1
        A(:,i) = 1;
    else
        A(:,i) = sin(i.*acos(grid))./sin(acos(grid));
    end
end
end