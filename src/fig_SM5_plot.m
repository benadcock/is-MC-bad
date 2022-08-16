%--- Description ---%
%
% Filename: fig_SM5_plot.m
% Authors: Ben Adcock and Simone Brugiapaglia 
% Part of the paper "Is Monte Carlo a bad sampling strategy for learning
% smooth functions in high dimensions?"
%
% Description: Plots Figure SM5

clear all; close all; clc;
addpath(genpath('../utils'))

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

poly_list = {'chebyshev','legendre','chebyshev2nd'}; % use Chebyshev (1st kind), Legendre or Chebyshev (2nd kind) polynomials
d_vals = [1 8 16]; % values of d to use
m_max_to_use = 3000; % maximum size of m to use
num_trials = 50; % number of trials
K = 100000; % error grid size
scale_type = 'log'; % use the logarithmic scaling (4.1) 
samp_type = 'MC'; % use Monte Carlo sampling
fun_list = {'iso_exp','reciprocal_linear'}; % functions to approximate

row_num = 1;

for fun_name = fun_list
    
    col_num = 1;
    
    for poly_type = poly_list
        
        fig = figure(1);
        h = []; cs = 1;
        
        for d = d_vals
            
            file_name = ['ALS_',samp_type,'_',char(poly_type),'_',char(fun_name),'_d',num2str(d),'_scaling',scale_type,'_trials',num2str(num_trials),'_K',num2str(K)];
            load(['../data/',file_name]);
            
            X = n_vals_data;
            Y = kappa_data;
            
            n_vals_all = [];
            mean_vals_all = [];
            curve_min_vals_all = [];
            curve_max_vals_all = [];
            
            for n = 1:m_max
                
                % find the trials which hit this value of m
                I = find(X == n);
                data_n = Y(I);
                
                if isempty(data_n) == 0
                    
                    % compute geometric mean and standard deviations
                    geo_mean = 10^(mean(log10(data_n)));
                    geo_std = std(log10(data_n));
                    curve_min = 10^(log10(geo_mean) - geo_std);
                    curve_max = 10^(log10(geo_mean) + geo_std);
                    
                    mean_vals_all = [mean_vals_all geo_mean];
                    curve_min_vals_all = [curve_min_vals_all curve_min];
                    curve_max_vals_all = [curve_max_vals_all curve_max];
                    n_vals_all = [n_vals_all n ];
                    
                end
                
            end
            
            % plot curves bounding the shaded region
            plot(n_vals_all,curve_min_vals_all, 'color', [colors{cs}, AlphaLevel]);
            hold on
            plot(n_vals_all,curve_max_vals_all, 'color', [colors{cs}, AlphaLevel])
            
            % fill in shaded region
            n_vals_2 = [n_vals_all, fliplr(n_vals_all)];
            inBetween = [curve_min_vals_all, fliplr(curve_max_vals_all)];
            fill(n_vals_2, inBetween, colors{cs}, 'FaceAlpha', AlphaLevel, 'EdgeAlpha', 0);
            
            % plot mean curve
            h = [h plot(n_vals_all,mean_vals_all,'LineWidth',lw,'Color',colors{cs})];
            
            cs = cs+1;
            
        end
        
        if isequal(char(poly_type),'legendre')
            pow = 2;
        elseif isequal(char(poly_type),'chebyshev')
            pow = log(3)/log(2);
        else
            pow = 3;
        end
        
        h = [h plot(n_vals_all,n_vals_all,'--','LineWidth',lw,'Color','k')];
        h = [h plot(n_vals_all,n_vals_all.^pow,'-.','LineWidth',lw,'Color','k')];
        
        set(gca, 'xscale', 'log')
        set(gca, 'yscale', 'log')
        hold off
        axis tight
        
        if isequal(char(poly_type),'legendre')
            leg = legend(h,'$d = 1$','$d = 8$','$d = 16$','$n$','$n^2$','location','northwest');
        elseif isequal(char(poly_type),'chebyshev')
            leg = legend(h,'$d = 1$','$d = 8$','$d = 16$','$n$','$n^{\log(3)/\log(2)}$','location','northwest');
        else
            leg = legend(h,'$d = 1$','$d = 8$','$d = 16$','$n$','$n^3$','location','northwest');
        end
        set(leg,'Interpreter','latex');
        
        xlabel('$n$','interpreter','latex');
        ylabel('$\kappa(\mathcal{P}_S)$','interpreter','latex');
        
        set_axis_param
        set_fonts
        
        saveas(fig,['../figs/fig_SM5_',num2str(row_num),num2str(col_num)],'epsc');
        col_num = col_num + 1;
        
    end
    
    row_num = row_num + 1;
end