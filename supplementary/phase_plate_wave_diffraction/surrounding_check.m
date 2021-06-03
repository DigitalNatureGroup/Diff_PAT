clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates supplementary figure 6 in the manuscript from the dataset.

N = 256;
f0 = 2e06;
c0 = 1480;
lambda = c0/f0;
dx = 0.00015; %smallest resolution printable on formalbs 2  https://formlabs.com/blog/3d-printer-resolution-meaning/
z_pos = 0.02;

image_read = (imread('test_transducer.png'));
image_read = reshape(double(image_read(:)), [N, N]);
input_amp = (image_read ./ 255);

target_amp = 1;
ii = 1;
phase_image_dp = double(imread(['test_image_' num2str(ii) '_phase_ex_DP.png']))./255.*(2*pi);
phase_image_gs = double(imread(['test_image_' num2str(ii) '_phase_ex_IASA.png']))./255.*(2*pi);

dp_image_1 = input_amp.*exp(1j.*phase_image_dp);
dp_image_1 = padarray(dp_image_1, [1024 1024]); % Added the zero padding to include regions outside the transducer regions
dp_image_1 = angularSpectrumCW(dp_image_1, dx, z_pos, f0, c0);

gs_image_1 = input_amp.*exp(1j.*phase_image_gs);
gs_image_1 = padarray(gs_image_1, [1024 1024]); % Added the zero padding to include regions outside the transducer regions
gs_image_1 = angularSpectrumCW(gs_image_1, dx, z_pos, f0, c0);

figure
pcolor(flipud(abs(dp_image_1)))
shading interp
title('Diff-PAT')
caxis([0 1])
axis equal
set(gca,'FontSize',24)
set(gca, 'Visible', 'off')
colorbar
exportgraphics(gca,['output_image\DP' num2str(ii) '.png'],'Resolution',300)

clf
pcolor(flipud(abs(gs_image_1)))
shading interp
title('IASA')
caxis([0 1])
axis equal
set(gca,'FontSize',24)
set(gca, 'Visible', 'off')
colorbar
exportgraphics(gca,['output_image\GS' num2str(ii) '.png'],'Resolution',300)
