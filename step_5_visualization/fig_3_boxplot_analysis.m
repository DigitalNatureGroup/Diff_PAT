clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

%% Main part begins
nm_array = [2]; % picking control number 
id_list = [3]; % choosing CES solver

%% Loading up dataset
file_name = ['result_14x14x1\result-others_sum-simple\ComparisonVarA_ID3_'];
T = readtable([file_name num2str(2) '.csv'],'Delimiter',',');
A = table2array(T(:, 4)); % picking target amplitude
res = table2array(T(:,6));

T1 = reshape(res, [2 1000]);
A1 = reshape(A, [2 1000]);

% Calculating interquartile range
Y = quantile(res,[0.25 0.75]);
inter_q = Y(2)-Y(1);


%% Within Q1 and Q3
q1_q3_c1 = T1(1,:);
index_c1 = find(q1_q3_c1 > Y(1) & q1_q3_c1 < Y(2));
q1_q3_c2 = T1(2,:);
index_c2 = find(q1_q3_c2 > Y(1) & q1_q3_c2 < Y(2));

%% Find out control geometries where both of them are in Q1 and Q3
intersect_c1_c2 = intersect(index_c1, index_c2); 

%% Outside of Q1 and Q3
index_c_out_c1 = find(q1_q3_c1 < Y(1) | q1_q3_c1 > Y(2));
index_c_out_c2 = find(q1_q3_c2 < Y(1) | q1_q3_c2 > Y(2));

%% Find out control geometries where both of them are outside of Q1 and Q3
intersect_out_c1_c2 = intersect(index_c_out_c1, index_c_out_c2);

%% Find out amplitude difference
A_min = min(A1(:, intersect_c1_c2)); %In Q1 - Q3
A_max = max(A1(:, intersect_c1_c2));
A_in = abs(A_max - A_min);

A_min = min(A1(:, intersect_out_c1_c2)); %Out of Q1 - Q3
A_max = max(A1(:, intersect_out_c1_c2));
A_out = abs(A_max - A_min);

disp(['N = 2; Mean Min-Max distance within Q1 and Q3:' num2str(mean(A_in))])
disp(['N = 2; Mean Min-Max distance outside Q1 and Q3:' num2str(mean(A_out))])

%% See overall delta in control point N =2 to 32
nm_array = [2 4 8 16 32];
for ii = 1:length(nm_array)
    file_name = ['result_14x14x1\result-others_sum-simple\ComparisonVarA_ID3_'];
    T = readtable([file_name num2str(nm_array(ii)) '.csv'],'Delimiter',',');
    A = table2array(T(:, 4)); % picking target amplitude
    A1 = reshape(A, [nm_array(ii) 1000]);
    
    A_min = min(A1);
    A_max = max(A1);
    
    A_diff = abs(A_max-A_min);
    
    disp(['Control Point N = ' num2str(nm_array(ii)) ' Mean Min-Max Distance Overall: ' num2str(mean(A_diff))])
end



