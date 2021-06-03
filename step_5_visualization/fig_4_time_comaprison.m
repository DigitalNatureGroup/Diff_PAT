clear all
close all
clc


%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp


%% This portion calculates the interquartile range of the error from fig. 3
nm_array = [2 4 8 16 32];
id_list = [2 3 4 5]; % for ES, CES, GS, Diff-PAT, respectively

iter_pick = 150;
store_data = zeros(length(id_list), length(nm_array), 3);
std_data = zeros(length(id_list), length(nm_array), 3);
iq_data = zeros(length(id_list), length(nm_array), 3);
for kk = 1:3
    for ii = 1:length(nm_array)
        for jj = 1:length(id_list)
            if kk == 1
                if id_list(jj) ~= 5
                    file_name = ['result_14x14x1\other_calc_time_lsit-96main\CalcTimeList_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    if nm_array(ii) < 10
                        file_name = ['result_14x14x1\our_calc_time_list\calc_time_list_00' num2str(nm_array(ii))];
                    else
                        file_name = ['result_14x14x1\our_calc_time_list\calc_time_list_0' num2str(nm_array(ii))];
                    end
                end
            elseif kk == 2
                if id_list(jj) ~= 5
                    file_name = ['result_16x16x2\other_calc_time_lsit-96main\CalcTimeList_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    if nm_array(ii) < 10
                        file_name = ['result_16x16x2\our_calc_time_list\calc_time_list_00' num2str(nm_array(ii))];
                    else
                        file_name = ['result_16x16x2\our_calc_time_list\calc_time_list_0' num2str(nm_array(ii))];
                    end
                end
            else
                if id_list(jj) ~= 5
                    file_name = ['result_32x32x1\other_calc_time_lsit-96main\CalcTimeList_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    if nm_array(ii) < 10
                        file_name = ['result_32x32x1\our_calc_time_list\calc_time_list_00' num2str(nm_array(ii))];
                    else
                        file_name = ['result_32x32x1\our_calc_time_list\calc_time_list_0' num2str(nm_array(ii))];
                    end
                end
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
                if id_list(jj) ~= 5
                    file_name = ['result_14x14x1\result-others_sum-simple\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    file_name = ['result_14x14x1\iter_phase_result\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '_iter' num2str(iter_pick)];
                end
            elseif kk == 2
                if id_list(jj) ~= 5
                    file_name = ['result_16x16x2\result-others_sum-simple\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    file_name = ['result_16x16x2\iter_phase_result\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '_iter' num2str(iter_pick)];
                end
            else
                if id_list(jj) ~= 5
                    file_name = ['result_32x32x1\result-others_sum-simple\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii))];
                else
                    file_name = ['result_32x32x1\iter_phase_result\ComparisonVarA_ID' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '_iter' num2str(iter_pick)];
                end
            end
            Ter = readtable([file_name '.csv'],'Delimiter',',');
            if id_list(jj) ~= 5
                Ter = table2array(Ter(:,end-1));
            else
                Ter = table2array(Ter(:,end));
            end
            eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) '=Ter;']);
            
            data_set = eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii))]);
            median_value = round(median(data_set), 3);
            Y = round(quantile(data_set,[0.25 0.75]), 3); % Finding Q1 and Q3
            IQ = Y(2)-Y(1); % Finding interquartile range
            
            iq_data(jj, ii, kk) = IQ;
        end
    end
end

%% Figure setup
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
hold on
line_prop = {'k-','k--','k-.'};
name_px = {'ES', 'CES', 'GS', 'DP'};
sub_label = {'a','b','c','d'};
color_spec = {'k','r','b'};
max_iq = 0.25;

%% Load up data for calculation time
for jj = 1:length(id_list)
    nexttile
    hold on
    title(name_px{jj})
    %% Draw execution time
    for kk = 1:3
        plot(nm_array, squeeze(store_data(jj, :, kk)),line_prop{kk},'LineWidth',2);
        set(gca, 'YScale', 'log')
    end
    
    %% Overlay interquartile range
    for kk = 1:3
        marker_col = (iq_data(jj, :, kk)./max_iq);
        scatter(nm_array, squeeze(store_data(jj, :, kk)), 30+marker_col*500, 1-marker_col,'filled', 'LineWidth',2,'MarkerEdgeColor','k');
        set(gca, 'YScale', 'log')
    end
    
    %% for illustrator file
    if jj == length(id_list)
        marker_col = linspace(0, 1, 5);
        scatter(nm_array, [10 10 10 10 10], 30+marker_col*500, 1-marker_col,'filled', 'LineWidth',2,'MarkerEdgeColor','k');
        colormap gray
    end
    
    %% Figure setup
    legend('M = 196','M = 512','M = 1024','Location','southeast')
    ylim([0.01 3000])
    set(gca,'FontSize',23)
    set(gca, 'FontName', 'Arial')
    if or(jj == 1, jj == 3)
        ylabel('Average Time [ms]')
    else
        set(gca,'yticklabel',[])
    end
    
    if or(jj ==3, jj == 4)
        xlabel('N [-]')
    else
        set(gca,'xticklabel',[])
    end
    colormap gray
    caxis([0 1])
    text(-5 ,10000, sub_label{jj},'FontName','Arial','FontSize',24,'FontWeight','Bold')
end
%% Rendering figures to PDF
fig= gcf;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 60 40],...
    'PaperSize',[60 40]);
print(fig,['time_graph'],'-dpdf','-r300')