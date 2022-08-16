%--- Description ---%
%
% Filename: generate_data_ALS_SM_scalings.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the data for the ALS scheme in the supplementaray
% materials for the experiment in Section SM4.2

clear all; close all; clc;
addpath(genpath('../utils'))

poly_type = 'legendre'; % use Legendre polynomials
d_vals = [1 2 4 8 16 32];
m_max = 1000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_list = {'linear1','linear2'}; % use the scalings (SM4.3)
samp_type = 'Opt'; % use the near-optimal sampling strategy
fun_list = {'iso_exp','sinusoid','reciprocal_linear','lognormal_ppde'}; % functions to approximate

for fun_name = fun_list
    for scale_type = scale_list
        for d = d_vals 
            run_ALS(poly_type,samp_type,char(fun_name),d,m_max,char(scale_type),num_trials,K) 
        end
    end
end