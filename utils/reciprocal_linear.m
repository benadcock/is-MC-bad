%--- Description ---%
%
% Filename: reciprocial_linear.m
% Authors: Ben Adcock and Simone Brugiapaglia
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the function f = f_3 defined in (4.4)
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = reciprocal_linear(y)

[m,d] = size(y);

b = ones(m,1);

for k = 1:d
    if d >= 2
        q = 10^(-3*(k-1)/(d-1));
    else
        q = 1;
    end
    b = b + (1/(2*d)).*y(:,k).*q;
end
b = 1./b;

end