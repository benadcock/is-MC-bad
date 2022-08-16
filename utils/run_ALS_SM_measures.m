%--- Description ---%
%
% Filename: run_ALS_SM_measures.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Runs the ALS approximation scheme for the different sampling measures
% and evaluates the error and condition number
%
% Inputs:
% poly_type - either 'legendre', 'chebyshev' or 'chebyshev2nd'
% samp_type - either 'MC' or 'Opt'
% fun_name - name of function, must have associated matlab function
% d - dimension
% m_max - maximum size of m to use
% scale_type - either 'log' (see (4.1)) or 'linear1' or 'linear2' (see
% (SM4.3))
% num_trials - number of trials
% K - error grid size

function run_ALS_SM_measures(poly_type,samp_type,fun_name,d,m_max,scale_type,num_trials,K)

%%% Define main parameters %%%

if isequal(scale_type,'log')
    scale_fun = @(t) max(t+1,ceil(t.*log(t))); % scaling to use between m and n
elseif isequal(scale_type,'linear1')
    scale_fun = @(t) ceil(1.5*t);
else
    scale_fun = @(t) ceil(2*t);
end

file_name = ['ALS_',samp_type,'_',poly_type,'_',fun_name,'_d',num2str(d),'_scaling',scale_type,'_trials',num2str(num_trials),'_K',num2str(K)];

space = ' ';

func = str2func(fun_name); % convert function name to function handle

beta = 0.5; % used in the bulk procedure

%%% Main loop %%%

err_grid = generate_sampling_grid(poly_type,d,K); % generate error grid
b_err_grid = func(err_grid)/sqrt(K); % generate function values over the error grid

% arrays for storing the data
err_data = zeros(num_trials,m_max+1);
cond_num_data = zeros(num_trials,m_max+1);
m_vals_data = zeros(num_trials,m_max+1);
n_vals_data = zeros(num_trials,m_max+1);
kappa_data = zeros(num_trials,m_max+1);

% loop over the trials
parfor t = 1:num_trials
    
    S = zeros(d,1);  % initialize the index set
    n = 1; u = 1;
    
    err_data_single = zeros(1,m_max+1);
    cond_num_data_single = zeros(1,m_max+1);
    m_vals_data_single = zeros(1,m_max+1);
    n_vals_data_single = zeros(1,m_max+1);
    kappa_data_single = zeros(1,m_max+1);
    
    % find the maximum n, given the scaling
    a = 1:m_max;
    a = scale_fun(a);
    n_max = max(find(a <= m_max));
    
    while n < n_max
        
        n = size(S,2);
        m = scale_fun(n); % define m
        
        n_vals_data_single(u) = n;
        m_vals_data_single(u) = m;
        
        A_err_grid = generate_measurement_matrix(poly_type,S,err_grid); % generate error matrix
        
        %%% generate measurement matrix, sample points and measurement vector %%%
        [Q,~] = qr(A_err_grid,0);
        
        if isequal(samp_type,'Opt')
            prob = sum(abs(Q).^2,2);
        else
            prob = ones(K,1);
        end
        prob = prob/sum(prob);
        
        L = datasample((1:K)',m,'Replace',true,'Weights',prob);
        W = diag(1./sqrt(K*prob(L)));
        y_grid = err_grid(L,:);
        
        %A = W*generate_measurement_matrix(poly_type,S,y_grid);
        A = sqrt(K/m)*W*Q(L,:);
        
        b = W*func(y_grid)/sqrt(m);
        
        %%% compute least-squares fit %%%
        c = A\b;
        
        %%%% evaluate on error grid %%%
        %bapprox_err_grid = A_err_grid*c;
        bapprox_err_grid = Q*c;
        
        L2_err = norm(bapprox_err_grid - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
        condA = cond(A);
        err_data_single(u) = L2_err;
        cond_num_data_single(u) = condA;
        
        
        disp(['ALS-',samp_type,': ',fun_name,space,' d = ',num2str(d),space,' trial = ',num2str(t),space,' n = ',num2str(n),space,' m = ',num2str(m),space,' Err = ',num2str(L2_err),' Cond = ',num2str(condA)]);
        
        if isequal(samp_type,'MC')
            w = generate_intrinsic_weights(poly_type,S);
            kappa_data_single(u) = sum(w.^2);
        end
        
        %%% Update the set S %%%
        
        RS = find_margin(S); % find the reduced margin of S
        
        % compute coefficients in reduced margin
        B = generate_measurement_matrix(poly_type,RS,y_grid);
        if isequal(samp_type,'MC')
            fcoeffs = (b-A*c)'*B;
        else
            fcoeffs = (b-A*c)'*W*B;
        end
        
        % bulk procedure
        fsum = sum(abs(fcoeffs).^2);
        [fcoeffs_sort,L] = sort(abs(fcoeffs).^2,'descend');
        
        fpsum = 0;
        i = 0;
        while fpsum < beta*fsum
            i = i+1;
            fpsum = fpsum + fcoeffs_sort(i);
        end
        
        F = RS(:,L(1:i));
        
        % update S and n
        S = [S F];
        n = size(S,2);
        
        u = u+1; % update counter
        
    end
    
    err_data(t,:) = err_data_single;
    cond_num_data(t,:) = cond_num_data_single;
    n_vals_data(t,:) = n_vals_data_single;
    m_vals_data(t,:) = m_vals_data_single;
    kappa_data(t,:) = kappa_data_single;
    
end

%%% Save data %%%
clear A B Q R err_grid A_err_grid b_err_grid
save(['../data/',file_name,'.mat'])

end


