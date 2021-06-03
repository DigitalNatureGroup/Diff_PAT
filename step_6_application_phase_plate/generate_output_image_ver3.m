clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates Fig. 4 in the main manuscript from the dataset.

N = 256; % Resolution
f0 = 2e06; % acoustic frequency Hz
c0 = 1480; % Speed of sound in water m/s
lambda = c0/f0; % Acoustic wavelength
dx = 0.00015; % Print Resolution
z_pos = 0.02; % XY plane position in Z axis

% Loading up transducer footprint
image_read = (imread('generating_inputs\test_package\test_transducer.png'));
image_read = reshape(double(image_read(:)), [N, N]);
input_amp = (image_read ./ 255);

number_of_samples = 5;
target_amp = 1;

for ii = 1:number_of_samples
    %%Load up target image
    test_image_original = double(imread(['generating_inputs\test_package\test_image_' num2str(ii) '_amp.png']));
    test_image_original = test_image_original ./ 255;
    test_image_original = test_image_original ./ max(max(test_image_original)) .* target_amp;
    
    %% Load up optimized acoustic holograms from each optimizer
    phase_image_dp = double(imread(['Diff_PAT\test_image_' num2str(ii) '_phase_ex_DP.png']))./255.*(2*pi);
    phase_image_gs = double(imread(['IASA\test_image_' num2str(ii) '_phase_ex_IASA.png']))./255.*(2*pi);
    
    % Use angualar spectrum to propagage the acoustic holgoram. Requires
    % k-Wave to do this. 
    dp_image_1 = input_amp.*exp(1j.*phase_image_dp);
    dp_image_1 = abs(angularSpectrumCW(dp_image_1, dx, z_pos, f0, c0));
    
    gs_image_1 = input_amp.*exp(1j.*phase_image_gs);
    gs_image_1 = abs(angularSpectrumCW(gs_image_1, dx, z_pos, f0, c0));
    
    % Use PSNR function from Image Processing Toolbox to evaluate the
    % peak-to-signal ratio
    dp_i = psnr(dp_image_1, test_image_original);
    gs_i = psnr(gs_image_1, test_image_original);
    
    % Render target image
    figure(1)
    clf
    pcolor(flipud(test_image_original))
    shading interp
    axis equal
    caxis([0 target_amp])
    set(gca, 'Visible', 'off')
    ax = gca;
    % Requires R2020a or later
    exportgraphics(ax,['output_image\original' num2str(ii) '.png'],'Resolution',300)
    
    % Render Diff-PAT image
    figure(1)
    clf
    pcolor(flipud(dp_image_1))  
    shading interp
    axis equal
    caxis([0 target_amp])
    set(gca, 'Visible', 'off')               
    ax = gca;
%     Requires R2020a or later
    exportgraphics(ax,['output_image\DP' num2str(ii) '.png'],'Resolution',300)

    % Render IASA image
    figure(1)
    clf
    pcolor(flipud(gs_image_1))
    hold on
    mx_gs = abs(test_image_original-gs_image_1);
    max_diff_gs = max(max(mx_gs));
    [row, column] = find(flipud(mx_gs) == (max_diff_gs));
    scatter(column, row, 120, 'rx','LineWidth',2);
    
    axis equal
    shading interp
    caxis([0 target_amp])
    set(gca, 'Visible', 'off')
    ax = gca;
    % Requires R2020a or later
    exportgraphics(ax,['output_image\IASA' num2str(ii) '.png'],'Resolution',300)
    
    disp(['Image ' num2str(ii) ' ---  DP: ' num2str(dp_i) ' GS: ' num2str(gs_i)])
    disp(['Max DP deviation: ' num2str(max_diff_dp)])
    disp(['Max GS deviation: ' num2str(max_diff_gs)])
    disp(['Max GS: ' num2str(max(max(abs(gs_image_1))))])
end