%--- Description ---%
%
% Filename: define_diffusion_term_1D.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: computes the diffusion coefficient a(x,y) of the parametric
% DE defined in Sec SM3.3
%
% Input:
% y - d x 1 sample point
%
% Output:
% b - m x 1 array of function values at the sample points

function a = define_diffusion_term_1D(y)

% Preliminary checks
if sum(y > 1) + sum(y < -1) > 0
    error('Parameter value is not in [-1,1]^d')
end
if min(size(y)) > 1
    error('y is not a vector')
end

% Extract parametric dimension
d = length(y);

% Define diffusion coefficient 
beta_c = 1/8;
beta_p = max([1,2*beta_c]);
beta = beta_c/beta_p; 
zeta = @(i) (sqrt(pi)*beta)^(1/2) * exp(-(floor(i/2)*pi*beta)^2/8);
a_exponent = @(x) 1 + y(1) .* (sqrt(pi).*beta/2).^(1/2);
for i = 2:d
    if mod(i,2)==0
        a_exponent = @(x) a_exponent(x) + zeta(i) .* sin(floor(i/2)*pi*x/beta_p) .* y(i);
    else
        a_exponent = @(x) a_exponent(x) + zeta(i) .* cos(floor(i/2)*pi*x/beta_p) .* y(i);
    end
end
a = @(x) exp(a_exponent(x));
