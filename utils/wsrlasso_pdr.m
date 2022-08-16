%--- Description ---%
%
% Filename: wsrlasso_pdr.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: solves the (weighted) SR-LASSO problem using restarted
% primal dual iteration (Algorithm SM2.3)
%
% Inputs:
% A - m x N measurement matrix
% b - m x 1 measurement vector
% w - N x 1 vector of weights (optional)
% lambda - SR-LASSO parameter
% sigma - PDI stepsize parameter
% tau - PDI stepsize parameter
% T - number of PDI iterations
% zeta - restart parameter
% r - restart parameter
% s - restart parameter
%
% Outputs:
% c - N x 1 array, the solution to the SR-LASSO problem
% u - number of restarts used (maximum is 100);
% c_diff - the final difference in iterates
%
% Note: if w is set as [], then the unweighted SR-LASSO problem is solved

function [c,u,c_diff] = wsrlasso_pdr(A,b,w,lambda,sigma,tau,T,zeta,r,s)

[m,N] = size(A);

if isequal(w,[])
    w = ones(N,1);
end

w = w(:); % the weights must form a column vector

if lambda >= 0
    
    % initialization
    c = zeros(N,1);
    epsl = norm(b);
    
    max_rsts = 100; % maximum number of restarts
    tol = 10*zeta;
    
    c_diff = 20*zeta;
    
    xi = zeros(m,1);
    p = zeros(N,1);
    q = zeros(m,1);
    csc = zeros(N,1); bsc = zeros(m,1);
    
    u = 0; % restart counter
    while c_diff > tol && u <= 100
        
        epsl = r*(epsl+zeta);
        a = s*epsl;
        
        % run PDI
        csc = c/a; bsc = b/a;
        xi = zeros(m,1);
        
        for n = 1:T
        p = csc - tau*A'*xi;
        csc_new = max(abs(p)-tau*lambda*w,0).*sign(p);
        q = xi+sigma*A*(2*csc_new-csc)-sigma*bsc;
        xi = min(1,1/norm(q))*q;
        csc = csc_new;
        end
        
        c_new = a*csc;
        c_diff = norm(c_new - c);
        c = c_new;
        
        u = u+1;
        
    end
    
else
    error('LAMBDA must be nonnegative')
end

end
