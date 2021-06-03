clear all
close all
clc

%%
% Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
% Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
% Corresponding email: tfushimi@slis.tsukuba.ac.jp

% This program generates supplementary figure 5 in the manuscript from the dataset.

N = 256;
f0 = 2e06;
c0 = 1480;
lambda = c0/f0;
dx = 0.00015; %smallest resolution printable on formalbs 2  https://formlabs.com/blog/3d-printer-resolution-meaning/
z_pos = 0.02;

image_read = (imread('generating_inputs\test_package\test_transducer.png'));
image_read = reshape(double(image_read(:)), [N, N]);
input_amp = (image_read ./ 255);

number_of_samples = 1;
sigma_levels = linspace(0, 0.3, 30);
store_psnr = zeros(1,length(sigma_levels));
target_amp = 1;

%Making results reproduceable by setting rand generator
rng('default');

for ii = 1:number_of_samples
    test_image_original = double(imread(['generating_inputs\test_package\test_image_' num2str(ii) '_amp.png']));
    test_image_original = test_image_original ./ 255;
    test_image_original = test_image_original ./ max(max(test_image_original)) .* target_amp;
    
    phase_image_dp = double(imread(['Diff_PAT\test_image_' num2str(ii) '_phase_ex_DP.png']))./255.*(2*pi);
    
    sz = size(phase_image_dp);
    for jj = 1:length(sigma_levels)
        %% Normally distributed noise with varying sigma level.
        noise_adder = normrnd(1, sigma_levels(jj), sz);
        
        dp_image_1 = input_amp.*exp(1j.*(phase_image_dp.*noise_adder));
        dp_image_1 = abs(angularSpectrumCW(dp_image_1, dx, z_pos, f0, c0));
        
        dp_i = psnr(dp_image_1, test_image_original);
        
        figure(1)
        clf
        pcolor(flipud(dp_image_1))
        shading interp
        axis equal
        caxis([0 1])
        colorbar
        set(gca, 'Visible', 'off')
        ax = gca;
        exportgraphics(ax,['output_image\dp_' num2str(ii) '_noise_' num2str(jj) '.png'],'Resolution',300)
        
        store_psnr(jj) = dp_i;
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(sigma_levels, store_psnr, '-o','LineWidth',2)
plot([min(sigma_levels.*100) max(sigma_levels)], [store_psnr(1)-3 store_psnr(1)-3],'r','LineWidth',2)
xlabel('\sigma of Normal Distribution (\mu = 1) [-]')
ylabel('PSNR [dB]')
set(gca,'FontSize',20)
legend({'PSNR','3 dB Drop'})
exportgraphics(gcf,'manufacturing_tolereance.png')