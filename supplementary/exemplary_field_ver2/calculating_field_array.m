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
side_length = 0.08; % m
center_z = 0.1;

lambda = c0/f0;
k = 2*pi*f0/c0;
a = 5e-03; % Transducer Radius

%% Drawing circle
circle_draw_phi = linspace(0, 2*pi);
circle_x = a*cos(circle_draw_phi);
circle_y = a*sin(circle_draw_phi);

%% Transducer Position
file_name = 'tpos_16x16';

T = readtable([file_name '.csv'],'Delimiter',',');
T = table2array(T);
trans = T;
%% Targets
target_index = 1;

T = readtable(['amplitudes_002.csv'],'Delimiter',',');
T = table2array(T);
amplitude_target = T(target_index, 2:3);

T = readtable(['geometries_002.csv'],'Delimiter',',');
T = table2array(T);
geometry_target = T(target_index, 2:end);
g_x = geometry_target(1:3:end);
g_y = geometry_target(2:3:end);
g_z = geometry_target(3:3:end);

figure('units','normalized','outerposition',[0 0 1 1])

solver_list = {'GS','Long2','Long3','ours'};

for sl = 1:length(solver_list)
    %%
    switch sl
        case 1
            T = readtable('GS-PAT_ID4_imag_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            imag_p = T(target_index,:);
            
            T = readtable('GS-PAT_ID4_real_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            real_p = T(target_index,:);
            
            combined = complex(real_p, imag_p);
            A_o = abs(combined);
            phase_o = angle(combined);
        case 2
            T = readtable('Long_ID2_imag_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            imag_p = T(target_index,:);
            
            T = readtable('Long_ID2_real_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            real_p = T(target_index,:);
            
            combined = complex(real_p, imag_p);
            A_o = abs(combined);
            phase_o = angle(combined);
        case 3
            T = readtable('Long_ID3_imag_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            imag_p = T(target_index,:);
            
            T = readtable('Long_ID3_real_2.csv','Delimiter',',');
            T = table2array(T(:,1:end-1));
            real_p = T(target_index,:);
            
            combined = complex(real_p, imag_p);
            A_o = abs(combined);
            phase_o = angle(combined);
        case 4
            T = readtable('phases_002.csv','Delimiter',',');
            T = table2array(T);
            phase_o = T(target_index,:);
            A_o = ones(length(phase_o),1);
    end
    A_o(A_o>1) = 1;
    while min(phase_o) < 0
        phase_o = phase_o + 2*pi;
    end
    
    phase_o = mod(phase_o, 2*pi);
    
    clf
    for tr = 1:length(trans)
        if trans(tr, 3) > 0
            subplot(2,2,3)
            hold on
            patch(circle_x-trans(tr, 1), circle_y-trans(tr, 2), A_o(tr).*[1 1 1])
        else
            subplot(2,2,1)
            hold on
            patch(circle_x-trans(tr, 1), circle_y-trans(tr, 2), A_o(tr).*[1 1 1])
        end
        
        if trans(tr, 3) > 0
            subplot(2,2,4)
            hold on
            patch(circle_x-trans(tr, 1), circle_y-trans(tr, 2), (phase_o(tr)./(2*pi)).*[1 1 1])            
        else
            subplot(2,2,2)
            hold on
            patch(circle_x-trans(tr, 1), circle_y-trans(tr, 2), (phase_o(tr)./(2*pi)).*[1 1 1])
        end
    end
    subplot(2,2,1)
    title('Amplitude z = 0 m')
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    set(gca,'FontSize',20)
    colorbar
    subplot(2,2,3)
    title('Amplitude z = 0.2355 m')
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    set(gca,'FontSize',20)
    colorbar
    
    subplot(2,2,2)
    title('Phase  z = 0 m')
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    set(gca,'FontSize',20)
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = linspace(0, 1, 3) ; %Create 8 ticks from zero to 1
    cbh.TickLabels = {'0', '\pi', '2\pi'};
    colormap gray
    
    subplot(2,2,4)
    title('Phase  z = 0.2355 m')
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    set(gca,'FontSize',20)
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = linspace(0, 1, 3) ; %Create 8 ticks from zero to 1
    cbh.TickLabels = {'0', '\pi', '2\pi'};
    colormap gray
    
    exportgraphics(gcf,['export\array_' solver_list{sl} '.png'], 'Resolution',300);
    clf
    for ii = 1:2
        x_c = -side_length:lambda/30:side_length;
        y_c = g_y(ii);
        z_c = center_z-side_length:lambda/30:center_z+side_length;
        
        [XX, YY, ZZ] = ndgrid(x_c, y_c, z_c);
        
        p1 = zeros(length(x_c), length(y_c), length(z_c));
        
        parfor tr = 1:length(trans)
            trans_x = trans(tr, :);
            trans_q = [0 0 1];
            
            if trans(tr, 3) > center_z
                trans_q(3) = -1;
            end
            
            r_prep_x = XX-trans_x(1);
            r_prep_y = YY-trans_x(2);
            r_prep_z = ZZ-trans_x(3);
            
            R = sqrt((r_prep_x).^2 + (r_prep_y).^2 + (r_prep_z).^2);
            dotproduct = r_prep_x.*trans_q(1) + r_prep_y.*trans_q(2) + r_prep_z.*trans_q(3);
            theta = acos(dotproduct./R./sqrt(trans_q(1).^2+trans_q(2).^2+trans_q(3).^2));
            D = directivity_fun(k, a, theta);
            
            p1 = p1 + (p_0./R).* A_o(tr) .* D .* exp(1j.*(k.*R+phase_o(tr)));
        end
        subplot(1,2,ii)
        hold on
        pcolor(x_c, z_c, squeeze(abs(p1))')
        scatter(g_x(ii), g_z(ii), 300,'kx', 'LineWidth', 2)
        title(['Target Amp = ' num2str(amplitude_target(ii)) ' Pa, Y = ' num2str(g_y(ii))])
        shading interp
        colormap hot
        caxis([0 2000])
        xlabel('X [m]')
        ylabel('Z [m]')
        axis equal
        set(gca,'FontSize',20)
        colorbar
    end
    exportgraphics(gcf,['export\' solver_list{sl} '.png'], 'Resolution',300);
end