%--- Description ---%
%
% Filename: run_CS.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Runs the CS approximation scheme and evaluates the error
%
% Inputs:
% poly_type - either 'legendre', 'chebyshev' or 'chebyshev2nd'
% samp_type - either 'MC' or 'Opt'
% fun_name - name of function, must have associated matlab function
% d - dimension
% m_max - maximum size of m to use
% m_vals_des - desired number of m values, these will be equally spaced
% between m = 1 and m = m_max
% num_trials - number of trials
% K - error grid size

function run_CS(poly_type,samp_type,fun_name,d,m_max,m_vals_des,num_trials,K)

%%% Define main parameters %%%

file_name = ['CS_',samp_type,'_',poly_type,'_',fun_name,'_d',num2str(d),'_trials',num2str(num_trials),'_K',num2str(K)];

space = ' ';

func = str2func(fun_name); % convert function name to function handle

N_des = 10*m_max; % desired maximum size of HC index set

zeta = 1e-15; % primal-dual restart parameter
r = exp(-1); % primal-dual restart parameter

% set the m values for CS
m_vals = round(linspace(1,m_max,m_vals_des));
num_m = length(m_vals);

err_grid = generate_sampling_grid(poly_type,d,K); % precompute error grid

n_CS = find_order('HC',d,N_des); % find the largest HC index set of size <= 6*m_max
I_CS = generate_index_set('HC',d,n_CS); % generate the HC index set
N = size(I_CS,2);

if isequal(samp_type,'Opt')
    % pre-generate error grid, then form QR and clear
    A_err_grid = generate_measurement_matrix(poly_type,I_CS,err_grid); % generate CS error matrix
    [Q,~] = qr(A_err_grid,0);
    clear A_err_grid;
    %prob = max(abs(Q).^2,2); % construct optimal sampling distribution
    prob = sum(abs(Q).^2,2);
else
    prob = ones(K,1);
end
prob = prob/sum(prob);

weights = generate_intrinsic_weights(poly_type,I_CS); % generate the intrinsic weights

coeff_data = zeros(num_trials,num_m,N);
rsts_data = zeros(num_trials,num_m);
c_diff_data = zeros(num_trials,num_m);

%%% Main loop %%%

% loop over m values
for i = 1:num_m
    
    m = m_vals(i);
    coeff_data_single = zeros(num_trials,N);
    rsts_data_single = zeros(num_trials,1);
    c_diff_data_single = zeros(num_trials,1);
    
    
    % loop over the trials
    parfor t = 1:num_trials
        
        %%% generate measurement matrix, sample points and measurement vector %%%
        
        L = datasample((1:K)',m,'Replace',true,'Weights',prob);
        W = diag(1./sqrt(K*prob(L)));
        y_grid = err_grid(L,:);
        A = W*generate_measurement_matrix(poly_type,I_CS,y_grid);
        b = W*func(y_grid)/sqrt(m);
        
        %%% solve using restarted primal-dual iteration %%%
        lambda = 1/sqrt(25*m); % value of lambda for wsrlasso
        normA = norm(A);
        sigma = 1/normA;
        tau = 1/normA;
        T = ceil(4*normA/r);
        s = T/(2*normA);
        
        [c,rsts,c_diff] = wsrlasso_pdr(A,b,weights,lambda,sigma,tau,T,zeta,r,s);
        
        coeff_data_single(t,:) = c;
        rsts_data_single(t) = rsts;
        c_diff_data_single(t) = c_diff;
        
        disp(['CS-',samp_type,': ','function = ',fun_name,space,' d = ',num2str(d),space,' trial = ',num2str(t),space,' m = ',num2str(m),space,' restarts  = ',num2str(rsts),' c_diff = ',num2str(c_diff)]);
        
    end
    
    coeff_data(:,i,:) = coeff_data_single;
    rsts_data(:,i) = rsts_data_single;
    c_diff_data(:,i) = c_diff_data_single;
    
end

%%% Generate error grid and evaluate errors %%%

A_err_grid = generate_measurement_matrix(poly_type,I_CS,err_grid); % generate CS error matrix
b_err_grid = func(err_grid)/sqrt(K);

err_data = zeros(num_trials,num_m);

for i = 1:num_m
    for t = 1:num_trials
        c = coeff_data(t,i,:); c = c(:);
        L2_err = norm(A_err_grid*c-b_err_grid)/norm(b_err_grid);
        err_data(t,i) = L2_err;
    end
end

%%% Save data %%%
clear A B Q R err_grid A_err_grid b_err_grid coeff_data coeff_data_single
save(['../data/',file_name,'.mat'])

end




