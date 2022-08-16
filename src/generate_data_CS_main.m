%--- Description ---%
%
% Filename: generate_data_CS_main.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Generates the data for the CS scheme in the main paper

clear all; close all; clc;
addpath(genpath('../utils'))

poly_type = 'legendre'; % use Legendre polynomials
d_vals = [1 2 4 8 16 32]; % values of d to use
m_max = 1000; % maximum size of m to use
m_vals_des = 20; % desired number of m values to use
num_trials = 50; % number of trials
K = 100000; % error grid size
samp_list = {'MC','Opt'}; % use Monte Carlo or the near-optimal sampling strategy
fun_list = {'iso_exp','sinusoid','reciprocal_linear'}; % functions to approximate

for fun_name = fun_list
    for samp_type = samp_list
        for d = d_vals
            run_CS(poly_type,char(samp_type),char(fun_name),d,m_max,m_vals_des,num_trials,K); 
        end
    end
end