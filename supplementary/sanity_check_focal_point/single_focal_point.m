clear all
close all
clc
%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program calculates the thoeretical maximum amplitude using
% analytical solution. Then outputs it as a target amplitude and
% geometries. Then use sanity_check_amplitude.py to evaluate the 80, 90,
% 110, 120 % of this theoretical maximum amplitude. 

p_0 = 1.98; % 12 V Pk-Pk Hirayama et al.
f0 = 40e03; % Hz
c0 = 346; % m/s
side_length = 0.05; % m

lambda = c0/f0;
k = 2*pi*f0/c0;
a = 5e-03; % Transducer Radius

T = readtable(['tpos_16x16.csv'],'Delimiter',',');
T = table2array(T);

center_z = max(T(:,3)) / 2;

%% Target points
x = linspace(-0.05, 0.05, 250);
y = zeros(length(x), 1);
z = ones(length(x), 1).*(center_z);

target_pressure = zeros(1, length(z));
A = zeros(250, 4);

for sample = 1:length(z)
    A(sample, :) = [sample-1 x(sample) y(sample) z(sample)];
    
    T_pos_tr = [x(sample); y(sample); z(sample)];
    p1 = 0;
    
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
    target_pressure(sample) = abs(p1);
end

B = [0:249; target_pressure]';
writematrix(A,'geometries_001.csv')
writematrix(B,'amplitudes_001.csv')