%--- Description ---%
%
% Filename: lognormal_ppde.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the parametric DE function defined in Sec SM3.3
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = lognormal_ppde(y)

g_FEM = @(x) 10 * ones(size(x)); % define forcing term (FEM format)
n = 1024 - 1; % DOF for FEM
[m,d] = size(y);

b = zeros(m,1);
for i = 1:m
    z = (y(i,:))';
    a = define_diffusion_term_1D(z);
    b(i) = FEM_1D_diffusion_QoI(a,g_FEM,n);
end