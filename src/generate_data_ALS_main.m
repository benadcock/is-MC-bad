%--- Description ---%
%
% Filename: generate_data_ALS_main.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the data for the ALS scheme in the main paper

clear all; close all; clc;
addpath(genpath('../utils'))

poly_type = 'legendre'; % use Legendre polynomials
d_vals = [1 2 4 8 16 32]; % values of d to use
m_max = 3000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 
samp_list = {'MC','Opt'}; % use Monte Carlo or the near-optimal sampling strategy
fun_list = {'iso_exp','sinusoid','reciprocal_linear','oned_fun'}; % functions to approximate

for fun_name = fun_list
    for samp_type = samp_list
        for d = d_vals 
            run_ALS(poly_type,char(samp_type),char(fun_name),d,m_max,scale_type,num_trials,K) 
        end
    end
end