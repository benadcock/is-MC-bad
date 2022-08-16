%--- Description ---%
%
% Filename: iso_exp.m
% Authors: Ben Adcock and Simone Brugiapaglia
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the function f = f_1 defined in (4.2)
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = iso_exp(y)

[m,d] = size(y);
b = exp(-sum(y,2)/(2*d));

end