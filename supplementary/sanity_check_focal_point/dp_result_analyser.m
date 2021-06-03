clear all
close all
clc
%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates supplementary figure 1 in the manuscript from the dataset.

p_0 = 1.98; % 12 V Pk-Pk Hirayama et al.
f0 = 40e03; % Hz
c0 = 346; % m/s
side_length = 0.05; % m

lambda = c0/f0;
k = 2*pi*f0/c0;
a = 5e-03; % Transducer Radius

T = readtable(['tpos_16x16.csv'],'Delimiter',',');
transducer_pos = table2array(T);

center_z = max(transducer_pos(:,3)) / 2;

%% Target points
x = linspace(-0.05, 0.05, 250);
y = zeros(length(x), 1);
z = ones(length(x), 1).*(center_z);

T = readtable('amplitudes_001.csv');
amp = table2array(T);
amp = amp(:,end);

percentage_m =[0.8,0.9,1.0, 1.1,1.2];
calculated_pressure = zeros(length(x), 1);

figure
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

for jj = 1:length(percentage_m)
    
    %% Draw solid line for 100 % amplitude
    if percentage_m(jj) == 1
        subplot(2,1,1)
        plot(x.*1000, amp.*percentage_m(jj),'k-','LineWidth',2.5);
        subplot(2,1,2)
        plot(x.*1000, amp.*percentage_m(jj),'k-','LineWidth',2.5);
    else % Otherwise calculate the field using phases stored in dp_results
        subplot(2,1,1)
        plot(x.*1000, amp.*percentage_m(jj),'LineWidth',2);
        
        T = readtable(['dp_results\phases_' num2str(round(percentage_m(jj)*100)) '.csv']);
        phases = table2array(T);
        for sample = 1:length(z)
            T_pos_tr = [x(sample); y(sample); z(sample)];
            p1 = 0;
            parfor tr = 1:length(transducer_pos)
                trans_x = [transducer_pos(tr, 1) transducer_pos(tr, 2) transducer_pos(tr, 3)];
                trans_q = [0 0 1];
                
                if transducer_pos(tr, 3) > center_z
                    trans_q(3) = -1;
                end
                
                r_prep_x = T_pos_tr(1)-trans_x(1);
                r_prep_y = T_pos_tr(2)-trans_x(2);
                r_prep_z = T_pos_tr(3)-trans_x(3);
                
                R = sqrt((r_prep_x).^2 + (r_prep_y).^2 + (r_prep_z).^2);
                dotproduct = r_prep_x.*trans_q(1) + r_prep_y.*trans_q(2) + r_prep_z.*trans_q(3);
                theta = acos(dotproduct./R./sqrt(trans_q(1).^2+trans_q(2).^2+trans_q(3).^2));
                D = directivity_fun(k, a, theta);
                
                p1 = p1 + (p_0./R) .* D .* exp(1j.*(k.*R+phases(sample,tr)));
            end
            calculated_pressure(sample) = abs(p1);
        end
        subplot(2,1,2)
        plot(x.*1000, calculated_pressure,'LineWidth',2);
    end
end
%% Rendering figure
subplot(2,1,1)

title('Target Amplitude, Y =0 m, Z = 0.1178 m')
grid on
grid minor
set(gca,'FontSize',14)
ylabel('Pressure [Pa]')
xlabel('Control Point Position X [mm]')
legend('80 %', '90 %', '100 %','110 %', '120 %','Location','eastoutside')
ylim([3500 6500])
subplot(2,1,2)
title('Optimized Solution using Diff-PAT')
legend('80 %', '90 %', '100 %','110 %', '120 %','Location','eastoutside')
set(gca,'FontSize',14)
ylabel('Pressure [Pa]')
xlabel('Control Point Position X [mm]')
grid on
grid minor
ylim([3500 6500])

exportgraphics(gcf,'analytical_versus_optimized.png')