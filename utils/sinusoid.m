%--- Description ---%
%
% Filename: sinusoid.m
% Authors: Ben Adcock and Simone Brugiapaglia
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the function f = f_2 defined in (4.3)
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = sinusoid(y)

[m,d] = size(y);

b = zeros(m,1);

for k = 1:d
    b = b + 0.3 + sin(16*y(:,k)/15 - 0.7) + (sin(16*y(:,k)/15 - 0.7)).^2 ;
end

end