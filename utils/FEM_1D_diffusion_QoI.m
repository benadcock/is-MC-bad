%--- Description ---%
%
% Filename: FEM_1D_diffusion_QoI.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: evaluates the QoI defined in Sec SM3.3
%
% Inputs:
% a - diffusion coefficient
% g - forcing term
% n - FEM DOF
%
% Output:
% QoI - the value of the QoI f = u(0.5)

function QoI = FEM_1D_diffusion_QoI(a,g,n)

h = 1/(n+1);

x_grid = @(i) i * h;

Q = zeros(6,n+1);
for i = 1:(n-1)
    Q(1,i) = 0;
    Q(2,i) = 0;
    Q(3,i) = 0;
    Q(4,i) = 1/(2 * h) * (    a(x_grid(i)) + a(x_grid(i-1)));
    Q(5,i) = 1/6 * h * (2 * g(x_grid(i)) + g(x_grid(i-1)));
    Q(6,i) = 1/6 * h * (2 * g(x_grid(i)) + g(x_grid(i+1)));
end
Q(2,n)   = 0;
Q(3,n)   = 0;
Q(4,n)   = 1/(2 * h) * (    a(x_grid(n))   + a(x_grid(n-1)));
Q(4,n+1) = 1/(2 * h) * (    a(x_grid(n+1)) + a(x_grid(n))  );
Q(5,n)   = 1/6 * h * (2 * g(x_grid(n))   + g(x_grid(n-1)));
Q(6,n)   = 1/6 * h * (2 * g(x_grid(n))   + g(x_grid(n+1)));


% Define stiffnedd matrix and loading vector
alpha = zeros(n,1);
beta = zeros(n-1,1);
b = zeros(n,1);
for i = 1:(n-1)
    alpha(i) = Q(4,i) + Q(4,i+1) + Q(2,i) + Q(3,i);
    beta(i)  = Q(1,i) - Q(4,i+1);
    b(i)     = Q(5,i) + Q(6,i);
end
alpha(n) = Q(4,n) + Q(4,n+1) + Q(2,n) + Q(3,n);
b(n)     = Q(5,n) + Q(6,n);

% Tridiagonal matrix
A = diag(beta,1);
A = A + A' + diag(alpha);

c = A \ b;

% QoI: approximation of u(0.5)
QoI = c(floor((n+1)/2));