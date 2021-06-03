clear all
close all
clc

%% 
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates Fig. 2 in the main manuscript from the dataset. 
% Tatsuki Fushimi

nm_array = [2 32]; % number of contrl points
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(length(nm_array),1, 'TileSpacing','Compact','Padding','Compact');
sampled_iter = [0:5:150]; % The iteration number plotted on the graph.
label_ix = {'a','b'};
%% Loop
for jj = 1:length(nm_array)
    nexttile
    hold on
    %%  Load up data for 14 x 14 array
    folder_name = ['result_14x14x1\iter_phase_result\ComparisonVarA_ID5_' num2str(nm_array(jj)) '_iter'];
    data_store = zeros(1000*nm_array(jj) , length(sampled_iter));
    for ii = 1:length(sampled_iter)
        T = readtable([folder_name num2str(sampled_iter(ii)) '.csv'],'Delimiter',',');
        T = table2array(T);
        data_store(:,ii) = T(:, end);
    end
    % Claculate mean and standard deviation
    mean_calc = mean(data_store);
    sd_calc = std(data_store);
    
    % Plot mean and standard deviation. S.D. as errorbar
    k1 = scatter(sampled_iter, mean_calc, 40, 'ko','filled','LineWidth',2);
    errorbar(sampled_iter, mean_calc, sd_calc,'k-','LineWidth', 2);
    ylim([0 2])
    
    %%  Load up data for 16 x 16 array
    folder_name = ['result_16x16x2\iter_phase_result\ComparisonVarA_ID5_' num2str(nm_array(jj)) '_iter'];
    data_store = zeros(1000*nm_array(jj) , length(sampled_iter));
    for ii = 1:length(sampled_iter)
        T = readtable([folder_name num2str(sampled_iter(ii)) '.csv'],'Delimiter',',');
        T = table2array(T);
        data_store(:,ii) = T(:, end);
    end
    
    % Claculate mean and standard deviation
    mean_calc = mean(data_store);
    sd_calc = std(data_store);
    
    % Plot mean and standard deviation. S.D. as errorbar
    k2 = scatter(sampled_iter, mean_calc, 40, 'ro','filled','LineWidth',2);
    errorbar(sampled_iter, mean_calc, sd_calc,'r-','LineWidth', 2);
    ylim([0 2])    
    
    %%  Load up data for 32 x 32 array
    folder_name = ['result_32x32x1\iter_phase_result\ComparisonVarA_ID5_' num2str(nm_array(jj)) '_iter'];
    data_store = zeros(1000*nm_array(jj) , length(sampled_iter));
    for ii = 1:length(sampled_iter)
        T = readtable([folder_name num2str(sampled_iter(ii)) '.csv'],'Delimiter',',');
        T = table2array(T);
        data_store(:,ii) = T(:, end);
    end
    
    % Claculate mean and standard deviation
    mean_calc = mean(data_store);
    sd_calc = std(data_store);
    
    % Plot mean and standard deviation. S.D. as errorbar
    k3 = scatter(sampled_iter, mean_calc, 40, 'bo','filled','LineWidth',2);
    errorbar(sampled_iter, mean_calc, sd_calc,'b-','LineWidth', 2);
    ylim([0 2])    
    %% Labeling and graphics touches
    
    if jj ~= length(nm_array)
        set(gca,'xticklabel',[])
    else
        xlabel('Iteration #')
    end
    ylabel('R_P [-]', 'FontSize', 24)
    set(gca,'FontSize',24)
    title(['N = ' num2str(nm_array(jj))])
    
    legend([k1, k2, k3], {'M = 196', 'M = 512', 'M = 1024'})
%     text(-18, 2.2, label_ix{jj}, 'FontName', 'Arial', 'FontSize', 24,'FontWeight', 'bold')
    set(gca, 'FontName', 'Arial')
end

%% Render figure as PDF
fig= gcf;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 30 20],...
    'PaperSize',[30 20]);
print(fig,['convergence_plot'],'-dpdf','-r0')