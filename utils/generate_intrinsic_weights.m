%--- Description ---%
%
% Filename: generate_intrinsic_weights.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: generates the intrinsic weights for polynomial approximation
%
% Inputs:
% poly_type - either 'legendre' (Legendre polynomials), 'chebyshev'
% (Chebyshev polynomials) or 'chebyshev2nd' (second kind Chebyshev polynomials)
% I - d x N array of multi-indices
%
% Output:
% u - N x 1 array of weights

function u = generate_intrinsic_weights(poly_type,I)

[d,N] = size(I); % get N (number of matrix columns) and d (dimension)
u = zeros(1,N);

if isequal(poly_type,'legendre')
    
    if d == 1
        u = sqrt(2.*I+ones(size(I)));      
    else
        u = prod(sqrt(2.*I+ones(size(I))));
    end
    
elseif isequal(poly_type,'chebyshev')
    
    for j = 1:N
        u(1,j) = sqrt(2)^(nnz(I(:,j)));
    end
    
elseif isequal(poly_type,'chebyshev2nd')
    
    if d == 1
        u = I+ones(size(I));      
    else
        u = prod(I+ones(size(I)));
    end
        
else
    error('Invalid poly_type');
end

u = u';
end


