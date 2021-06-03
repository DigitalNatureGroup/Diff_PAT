clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates Fig. 2 in the main manuscript from the dataset.

%% Picking number of control points and algorithm
nm_array = [2 4 8 16 32];
id_list = [2 3 4 5]; % for ES, CES, GS, Diff-PAT, respectively

% For graphics

x_pos_text = [1.2 2.2 3.2 4.2];
y_pos_text = 0.3;
x_place = -0.3;

% Which iteration number to pick for Diff-PAT
iter_pick = 150;
colors = fliplr({'b','y','r','g','k','b','y','r','g','k','b','y','r','g'});
figure('units','normalized','outerposition',[0 0 1 1])
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
% for control numbers
for ii = 1:length(nm_array)
    subplot(5,1, ii);
    hold on
    g = [];
    for kk = 1:3
        for jj = 1:length(id_list)
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
            T = readtable([file_name '.csv'],'Delimiter',',');
            if id_list(jj) ~= 5
                T = table2array(T(:,end-1));
            else
                T = table2array(T(:,end));
            end
            eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk) '=T;']);
        end
        
        % This generates data for xlabel
        g2 = repmat({['ES' num2str(kk)]},length(eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk)])),1);  %ID2
        g3 = repmat({['CES' num2str(kk)]},length(eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk)])),1); %ID3
        g4 = repmat({['GS' num2str(kk)]},length(eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk)])),1);  %ID4
        g5 = repmat({['DP' num2str(kk)]},length(eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk)])),1);  %ID5
        gspace = repmat({['  ' num2str(kk)]},length(eval(['id' num2str(id_list(jj)) '_' num2str(nm_array(ii)) 'M_' num2str(kk)])),1);  %ID5
        if kk ~= 3
            g = [g; g2; g3; g4; g5; gspace];
        else
            g = [g; g2; g3; g4; g5];
        end
    end
    % Plot box plot
    m = boxplot(eval(['[id2_' num2str(nm_array(ii)) 'M_1' ', id3_' num2str(nm_array(ii)) ...
        'M_1' ', id4_' num2str(nm_array(ii)) 'M_1' ', id5_' num2str(nm_array(ii)) 'M_1, -100.*ones(nm_array(ii)*1000, 1),'...
        'id2_' num2str(nm_array(ii)) 'M_2' ', id3_' num2str(nm_array(ii)) ...
        'M_2' ', id4_' num2str(nm_array(ii)) 'M_2' ', id5_' num2str(nm_array(ii)) 'M_2, -100.*ones(nm_array(ii)*1000, 1),'...
        'id2_' num2str(nm_array(ii)) 'M_3' ', id3_' num2str(nm_array(ii)) ...
        'M_3' ', id4_' num2str(nm_array(ii)) 'M_3' ', id5_' num2str(nm_array(ii)) 'M_3'...
        ']']),g,'Whisker',1.5,'Symbol','ko');
    hold on
    set(m,{'linew'},{2})
    plot([0 200],[1 1],'k:','LineWidth', 0.5) % Show where R_p = 1
    %% Graphics command
    ylim([0.65 1.2])
    ylabel('R_P [-]', 'FontSize', 24)
    
    if ii ~= length(nm_array)
        set(gca,'xticklabel','')
    else
        xticklabels({'ES','CES','GS','DP','','ES','CES','GS','DP','','ES','CES','GS','DP'})
    end
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', 34)
    h = findobj(gca,'Tag','Box');
    
    for kj=1:length(h)
        patch(get(h(kj),'XData'),get(h(kj),'YData'),colors{kj},'FaceAlpha',.5);
    end

end
%% Rendering figures to PDF
fig= gcf;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 60 60],...
    'PaperSize',[60 60]);
print(fig,['box_plots_exported\all'],'-dpdf','-r300')

