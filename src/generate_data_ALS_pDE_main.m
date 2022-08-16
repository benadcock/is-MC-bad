%--- Description ---%
%
% Filename: generate_data_ALS_pDE_main.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the data for the ALS scheme in the main paper for
% the parametric DE example

clear all; close all; clc;
addpath(genpath('../utils'))

poly_type = 'legendre'; % use Legendre polynomials
d_vals = [1 2 4 8 16 32]; % values of d to use
m_max_vals = [1000 1000 3000 1000 1000 3000]; % maximum sizes of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 
samp_list = {'MC','Opt'}; % use Monte Carlo or the near-optimal sampling strategy
fun_name = 'lognormal_ppde'; % function to approximate

for samp_type = samp_list
    for i = 1:6
        d = d_vals(i); m_max = m_max_vals(i);
        run_ALS(poly_type,char(samp_type),fun_name,d,m_max,scale_type,num_trials,K)
    end
end