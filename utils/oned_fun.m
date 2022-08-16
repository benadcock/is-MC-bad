%--- Description ---%
%
% Filename: oned_fun.m
% Authors: Ben Adcock and Simone Brugiapaglia
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the function f considered in Fig. 6
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = oned_fun(y)

[m,d] = size(y);

b = 1./(10-9*y(:,1));

end
