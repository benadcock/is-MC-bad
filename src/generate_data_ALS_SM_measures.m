%--- Description ---%
%
% Filename: generate_data_ALS_SM_measures.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the data for the ALS scheme in the supplementaray
% materials for the experiment in Section SM4.1

clear all; close all; clc;
addpath(genpath('../utils'))

poly_list = {'chebyshev','chebyshev2nd'}; % use Chebyshev (1st kind) or Chebyshev (2nd kind) polynomials
d_vals = [1 8 16]; % values of d to use
m_max = 3000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 
samp_list = {'MC','Opt'}; % use Monte Carlo or the near-optimal sampling strategy
fun_list = {'iso_exp','reciprocal_linear'}; % functions to approximate

for fun_name = fun_list
    for poly_type = poly_list
        for samp_type = samp_list
            for d = d_vals
                run_ALS_SM_measures(char(poly_type),char(samp_type),char(fun_name),d,m_max,scale_type,num_trials,K)
            end
        end
    end
end