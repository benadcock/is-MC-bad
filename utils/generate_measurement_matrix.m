%--- Description ---%
%
% Filename: generate_measurement_matrix.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: generates a measurement matrix using either tensor Chebyshev
% or Legendre polynomials from an arbitrary multi-index set and collection
% of sample points
%
% Inputs:
% poly_type - either 'legendre' (Legendre polynomials), 'chebyshev' (Chebyshev polynomials) or 'chebyshev2nd' (second kind Chebyshev polynomials)
% I - d x N array of multi-indices
% y_grid - m x d array of sample points
%
% Output:
% A - normalized measurement matrix A

function A = generate_measurement_matrix(poly_type,I,y_grid)

[d,N] = size(I); % get N (number of matrix columns) and d (dimension)
m = size(y_grid,1); % get m (number of matrix rows)
A = zeros(m,N); % initialize A

n = max(I(:)); % find maximum polynomial degree

for i = 1:m
    y = y_grid(i,:); % select ith sample point
    
    % evaluate the 1D polynomials at the components of y
    if isequal(poly_type,'legendre')
        L = legmat(y',n+1);
    elseif isequal(poly_type,'chebyshev')
        L = chebmat(y',n+1);
    elseif isequal(poly_type,'chebyshev2nd')
        L = chebmat2nd(y',n+1);
    else
        error('invalid poly_type')
    end
    
    % evaluate the dD polynomials via tensor products
    for j = 1:N
        Lij = zeros(d,1);
        for k = 1:d
            Lij(k,1) = L(k,I(k,j)+1);
        end
        A(i,j) = prod(Lij);
    end
end

% normalize A
A = A/sqrt(m);

end
