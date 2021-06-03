clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates supplementary figure 2 in the manuscript from the dataset.

index = 0:10:150;
N = [1 3 5];
figure('units','normalized','outerposition',[0 0 1 1])

for ii = 1:length(N)
    number_c = 2^N(ii);
    %% Load up loss function performance
    for jj = length(index)
        data = csvread(['original\ComparisonVarA_ID5_' num2str(number_c) '_iter' num2str(index(jj)) '.csv']);
        original_std_store = std(data(:,end));
        original_mean_store = mean(data(:,end));
        data = csvread(['loss_abs\ComparisonVarA_ID5_' num2str(number_c) '_iter' num2str(index(jj)) '.csv']);
        abs_std_store = std(data(:,end));
        abs_mean_store = mean(data(:,end));
        data = csvread(['loss_huber\ComparisonVarA_ID5_' num2str(number_c) '_iter' num2str(index(jj)) '.csv']);
        huber_std_store = std(data(:,end));
        huber_mean_store = mean(data(:,end));
    end
    subplot(1,3,ii)
    hold on
    %% Plot
    bar([1,2,3], [original_mean_store, huber_mean_store, abs_mean_store])
    er = errorbar([1,2,3], [original_mean_store, huber_mean_store, abs_mean_store], [original_std_store, huber_std_store, abs_std_store],'LineWidth',3);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    
    title(['N = ' num2str(number_c)])
    ylabel('R_p [-]')
    set(gca, 'XTick', [1 2 3])
    set(gca, 'XTickLabel', {'SSE' 'Huber','SAE'})
    ylim([0 1.5])
    set(gca,'FontSize',20)
end
%% Export
exportgraphics(gcf,'loss_comparison.png','Resolution',300)
