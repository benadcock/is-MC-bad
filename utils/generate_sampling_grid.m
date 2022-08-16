%--- Description ---%
%
% Filename: generate_sampling_grid.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: generates a random sampling grid using either the uniform, Chebyshev (1st kind) or Chebyshev (2nd kind) measure
%
% Inputs: 
% samp_type - either 'uniform' or 'legendre' (uniform measure), 'chebyshev' (Chebyshev
% measure) or 'chebyshev2nd' (Chebyshev second kind measure)
% d - dimension
% m - number of sample points
% 
% Output:
% y_grid - the m x d array where the ith row is the ith sample point y_i

function y_grid = generate_sampling_grid(samp_type,d,m)

% uniform measure
if isequal(samp_type,'uniform') || isequal(samp_type,'legendre')
    
y_grid = 2*rand(m,d)-ones(m,d); 

% Chebyshev measure
elseif isequal(samp_type,'chebyshev')

    y_grid = rand(m,d);
    y_grid = cos(pi .* y_grid); 

    
% Chebyshev 2nd kind measure
elseif isequal(samp_type,'chebyshev2nd')
    
    u = rand(m,d);
    v = rand(m,d);
    
    y_grid = sqrt(u).*cos(pi*v);
    
else

error('invalid samp_type')

end