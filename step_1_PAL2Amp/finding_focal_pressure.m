% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program calculates the maximum acoustic pressure amplitude possible
% by the given transducer arrays, "tpos_NxN.csv" file. 

clear all
close all
clc

%% Setting up acoustic field
p_0 = 1.98; % 12 V Pk-Pk Hirayama et al.
f0 = 40e03; % Hz
c0 = 346; % m/s
side_length = 0.05; % m

lambda = c0/f0;
k = 2*pi*f0/c0;
a = 5e-03; % Transducer Radius

%% Main Loop
for ii = 1:3
    % Loading CSV FIles
    if ii == 1
        file_name = 'tpos_14x14';
    elseif ii == 2
        file_name = 'tpos_16x16';
    elseif ii == 3
        file_name = 'tpos_32x32x1';
    end
    
    T = readtable([file_name '.csv'],'Delimiter',',');
    T = table2array(T);
    
    % Change Position of center of array in z depending on the array
    if or(ii == 1, ii == 3)
        center_z = 0.10;
    else
        center_z = max(T(:,3)) / 2;
    end
    
    % Get the location of the ROI, dependent on the center of array position. 
    x_c = [-side_length side_length];
    y_c = [-side_length side_length];
    z_c = [center_z-side_length center_z+side_length];
    
    [XX, YY, ZZ] = ndgrid(x_c, y_c, z_c);
    XX = XX(:);
    YY = YY(:);
    ZZ = ZZ(:);
    
    figure
    scatter3(XX, YY, ZZ)
    axis equal
    
    focal_pressure = zeros(1, length(XX));
    %% Loop to find acoustic pressure amplitude. Huygens' Principle
    % For every ROI point
    for sample = 1:length(XX)
        
        T_pos_tr = [XX(sample); YY(sample); ZZ(sample)];
        p1 = 0;
        % Go through the list of transducer, and linearly add them together
        for tr = 1:length(T)
            trans_x = [T(tr, 1) T(tr, 2) T(tr, 3)];
            trans_q = [0 0 1];
            
            if T(tr, 3) > center_z
                trans_q(3) = -1;
            end
            
            r_prep_x = T_pos_tr(1)-trans_x(1);
            r_prep_y = T_pos_tr(2)-trans_x(2);
            r_prep_z = T_pos_tr(3)-trans_x(3);
            
            R = sqrt((r_prep_x).^2 + (r_prep_y).^2 + (r_prep_z).^2);
            dotproduct = r_prep_x.*trans_q(1) + r_prep_y.*trans_q(2) + r_prep_z.*trans_q(3);
            theta = acos(dotproduct./R./sqrt(trans_q(1).^2+trans_q(2).^2+trans_q(3).^2));
            D = directivity_fun(k, a, theta);
            
            focal_phi = -(2*pi*f0/c0).*(sqrt((trans_x(1)-T_pos_tr(1)).^2+(trans_x(2)-T_pos_tr(2)).^2+(trans_x(3)-T_pos_tr(3)).^2) ...
                - sqrt(T_pos_tr(1).^2 + T_pos_tr(2).^2 + T_pos_tr(3).^2));
            p1 = p1 + (p_0./R) .* D .* exp(1j.*(k.*R+focal_phi));
        end
        focal_pressure(sample) = abs(p1);
    end
    %% Output
    format long
    disp(['Center of ROI is ' num2str(center_z) ' m'])
    disp(['Average Pressure at the edge of ROI in ' file_name ' is: ' num2str(round(mean(focal_pressure)))]);
end
