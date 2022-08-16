%--- Description ---%
%
% Filename: fig_SM1234_plot.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Plots Figures SM1-4

clear all; close all; clc;
addpath(genpath('../utils'))

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

poly_list = {'chebyshev','legendre','chebyshev2nd'}; % use Chebyshev (1st kind), Legendre or Chebyshev (2nd kind) polynomials
d_vals = [1 8 16]; % values of d to use
m_max_to_use = 3000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 
samp_list = fliplr({'MC','Opt'}); % use Monte Carlo or the near-optimal sampling strategy

fig_nums = [1 2 3 4];

for fig_num = fig_nums
    
    if fig_num == 1 || fig_num == 2
        fun_name = 'iso_exp';
    else
        fun_name = 'reciprocal_linear';
    end
    
    for poly_type = poly_list
        
        for d = d_vals
            
            fig = figure(1);
            h = []; cs = 1;
            
            for samp_type = samp_list
                
                file_name = ['ALS_',char(samp_type),'_',char(poly_type),'_',fun_name,'_d',num2str(d),'_scaling',scale_type,'_trials',num2str(num_trials),'_K',num2str(K)];
                load(['../data/',file_name]);
                
                X = m_vals_data;
                
                if fig_num == 1 || fig_num == 3
                    Y = err_data(:,:);
                else
                    Y = cond_num_data(:,:);
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
            
            if fig_num == 1 || fig_num == 3
                ylabel('Relative $L^2_{\varrho}(\mathcal{U})$-error','interpreter','latex');
            else
                ylabel('Condition number $\mathrm{cond}(\mathbf{A})$','interpreter','latex');
                if d >= 2 || isequal(poly_type,'chebyshev')
                    ylim([1e0,1e2]);
                else
                    ylim([1e0,1e17]);
                end
            end
            
            legend(h,'ALS-Opt','ALS-MC','location','northeast');
            
            set_axis_param
            set_fonts
            
            fig_name = ['fig_SM',num2str(fig_num),'_',num2str(poly_type),'_d',num2str(d)];
            saveas(fig,['../figs/',fig_name],'epsc');
            
        end
        
    end
end



