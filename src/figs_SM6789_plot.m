%--- Description ---%
%
% Filename: fig_SM6789_plot.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Plots Figures SM6-9

clear all; close all; clc;
addpath(genpath('../utils'))

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

poly_type = 'legendre'; % use Legendre polynomials
d_vals = [1 2 4 8 16 32]; % values of d to use
m_max_to_use = 1000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 

for fig_num = 6:9
    
    if fig_num == 6
        fun_name = 'iso_exp';
    elseif fig_num == 7
        fun_name = 'sinusoid';
    elseif fig_num == 8
        fun_name = 'reciprocal_linear';
    elseif fig_num == 9
        fun_name = 'lognormal_ppde';
    end
    
    for d = d_vals
        
        for fig_case = 1:2
            
            fig = figure(fig_case);
            h = [];
            
            if fig_case == 1
                csmax = 5;
            else
                csmax = 3;
            end
            
            for cs = 1:csmax
                
                if cs <= 3
                    
                    if cs == 1
                        scale_type = 'linear1';
                    elseif cs == 2
                        scale_type = 'linear2';
                    else
                        scale_type = 'log';
                    end
                    
                    file_name = ['ALS_Opt','_',poly_type,'_',fun_name,'_d',num2str(d),'_scaling',scale_type,'_trials',num2str(num_trials),'_K',num2str(K)];
                    load(['../data/',file_name]);
                    X = m_vals_data;
                    
                    if fig_case == 1
                        Y = err_data(:,:);
                    else
                        Y = cond_num_data(:,:);
                    end
                    
                    
                elseif cs == 4
                    file_name = ['CS_MC','_',poly_type,'_',fun_name,'_d',num2str(d),'_trials',num2str(num_trials),'_K',num2str(K)];
                    load(['../data/',file_name]);
                    X = m_vals; X = ones(num_trials,1).*m_vals;
                    Y = err_data;
                else
                    file_name = ['CS_Opt','_',poly_type,'_',fun_name,'_d',num2str(d),'_trials',num2str(num_trials),'_K',num2str(K)];
                    load(['../data/',file_name]);
                    X = m_vals; X = ones(num_trials,1).*m_vals;
                    Y = err_data;
                end
                
                m_vals_all = [];
                mean_vals_all = [];
                curve_min_vals_all = [];
                curve_max_vals_all = [];
                
                for m = 1:m_max_to_use
                    
                    % find the trials which hit this value of m
                    I = find(X == m);
                    data_m = Y(I);
                    
                    % compute the statistics of each data point
                    if isempty(data_m) == 0
                        
                        % compute geometric mean and standard deviations
                        geo_mean = 10^(mean(log10(data_m)));
                        geo_std = std(log10(data_m));
                        curve_min = 10^(log10(geo_mean) - geo_std);
                        curve_max = 10^(log10(geo_mean) + geo_std);
                        
                        mean_vals_all = [mean_vals_all geo_mean];
                        curve_min_vals_all = [curve_min_vals_all curve_min];
                        curve_max_vals_all = [curve_max_vals_all curve_max];
                        m_vals_all = [m_vals_all m ];
                        
                    end
                    
                end
                
                % plot curves bounding the shaded region
                plot(m_vals_all,curve_min_vals_all, 'color', [colors{cs}, AlphaLevel]);
                hold on
                plot(m_vals_all,curve_max_vals_all, 'color', [colors{cs}, AlphaLevel])
                
                % fill in shaded region
                m_vals_2 = [m_vals_all, fliplr(m_vals_all)];
                inBetween = [curve_min_vals_all, fliplr(curve_max_vals_all)];
                fill(m_vals_2, inBetween, colors{cs}, 'FaceAlpha', AlphaLevel, 'EdgeAlpha', 0);
                
                % plot mean curve
                h = [h plot(m_vals_all,mean_vals_all,'LineWidth',lw,'Color',colors{cs})];
                
                cs = cs+1;
            end
            
            set(gca, 'yscale', 'log')
            
            hold off
            axis tight
            
            xlabel('$m$','interpreter','latex');
            
            if fig_case == 1
                ylabel('Relative $L^2_{\varrho}(\mathcal{U})$-error','interpreter','latex');
                legend(h,'ALS-Opt, $1.5 n$','ALS-Opt, $2 n$','ALS-Opt, $n\log(n)$','CS-MC','CS-Opt','location','northeast');
            else
                ylabel('Condition number $\mathrm{cond}(\mathbf{A})$','interpreter','latex');
                legend(h,'ALS-Opt, $1.5 n$','ALS-Opt, $2 n$','ALS-Opt, $n\log(n)$','location','northeast');
                
                if d == 1
                    ylim([1,1e10]);
                else
                    ylim([1,1e2]);
                end
                
            end
            
            set_axis_param
            set_fonts
            
            if fig_case == 1
                fig_name = ['fig_SM',num2str(fig_num),'_d',num2str(d),'_err'];
            else
                fig_name = ['fig_SM',num2str(fig_num),'_d',num2str(d),'_cond'];
            end
            saveas(fig,['../figs/',fig_name],'epsc');
            
        end
        
    end
    
end


