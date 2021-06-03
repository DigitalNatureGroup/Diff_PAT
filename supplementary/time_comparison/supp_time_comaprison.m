clear all
close all
clc


%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates supplementary figure 3 in the manuscript from the dataset.


%% Picking number of control points and algorithm
nm_array = [2 4 8 16 32];
id_list = [5]; % for ES, CES, GS, Diff-PAT, respectively


iter_pick = 150;
store_data = zeros(length(id_list), length(nm_array), 2);
std_data = zeros(length(id_list), length(nm_array), 2);
iq_data = zeros(length(id_list), length(nm_array), 2);
for kk = 1:2
    for ii = 1:length(nm_array)
        for jj = 1:length(id_list)
            if kk == 1
                file_name = ['result_14x14x1\our_calc_time_list\calc_time_list_' num2str(nm_array(ii), '%03.f')];
            elseif kk == 2
                file_name = ['result_16x16x2\our_calc_time_list\calc_time_list_' num2str(nm_array(ii), '%03.f')];
            end
            T = readtable([file_name '.csv'],'Delimiter',',');
            T = table2array(T(:,end));
            if id_list(jj) == 5
                T = T*1000;
            end
            
            store_data(jj, ii, kk) = mean(T);
            std_data(jj, ii, kk) = std(T);
            %% calculating error
            if kk == 1
                file_name = ['result_14x14x1\iter_phase_result\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '_iter' num2str(iter_pick)];
                
            elseif kk == 2
                file_name = ['result_16x16x2\iter_phase_result\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '_iter' num2str(iter_pick)];
            end
            Ter = readtable([file_name '.csv'],'Delimiter',',');
            Ter = table2array(Ter(:,end));
            
            eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '=Ter;']);
            
            data_set = eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii))]);
            median_value = round(median(data_set), 3);
            Y = round(quantile(data_set,[0.25 0.75]), 3); % Finding Q1 and Q3
            IQ = Y(2)-Y(1); % Finding interquartile range
            
            iq_data(jj, ii, kk) = IQ;
        end
    end
end

julia_store_data = zeros(length(id_list), length(nm_array), 2);
julia_std_data = zeros(length(id_list), length(nm_array), 2);
julia_iq_data = zeros(length(id_list), length(nm_array), 2);
%% JULIA
for kk = 1:2
    for ii = 1:length(nm_array)
        for jj = 1:length(id_list)
            if kk == 1
                file_name = ['julia\14x14x1\timed_' num2str(nm_array(ii),'%03.f')];
            elseif kk == 2
                file_name = ['julia\16x16x2\timed_' num2str(nm_array(ii),'%03.f')];
            end
            
            T = readtable([file_name '.csv'],'Delimiter',',');
            T = table2array(T(:,1));
            if id_list(jj) == 5
                T = T*1000;
            end
            
            julia_store_data(jj, ii, kk) = mean(T);
            julia_std_data(jj, ii, kk) = std(T);
            %% calculating error
            if kk == 1
                file_name = ['julia\14x14x1\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
            elseif kk == 2
                file_name = ['julia\16x16x2\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
            end
            Ter = readtable([file_name '.csv'],'Delimiter',',');
            Ter = table2array(Ter(:,end));
            
            eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '=Ter;']);
            
            data_set = eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii))]);
            median_value = round(median(data_set), 3);
            Y = round(quantile(data_set,[0.25 0.75]), 3); % Finding Q1 and Q3
            IQ = Y(2)-Y(1); % Finding interquartile range
            
            julia_iq_data(jj, ii, kk) = IQ;
        end
    end
end


%% Plotting 
% Simply compare the computationl time. 
figure('units','normalized','outerposition',[0 0 1 1])
hold on
line_prop = {'k-+','k--+'};
name_px = {'DP'};

hold on

for kk = 1:2
    plot(nm_array, squeeze(store_data(jj, :, kk)),line_prop{kk},'LineWidth',2);
    set(gca, 'YScale', 'log')
end

line_prop = {'r-+','r--+'};
for kk = 1:2
    plot(nm_array, squeeze(julia_store_data(jj, :, kk)),line_prop{kk},'LineWidth',2);
    set(gca, 'YScale', 'log')
end

legend('M = 196 (AD)','M = 512 (AD)','M = 196 (ND)','M = 512 (ND)','Location','southeast')
set(gca,'FontSize',23)
set(gca, 'FontName', 'Arial')
ylabel('Average Time [ms]')
xlabel('N [-]')
colormap gray
caxis([0 1])

%% Rendering figures to PNG
exportgraphics(gcf, 'time_comparison_AD_ND.png','Resolution',300)