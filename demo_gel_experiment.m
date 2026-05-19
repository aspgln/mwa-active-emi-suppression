% =========================================================================
% Script: Demo_gel_experiment.m
% Author: Qing Dai
% Date: 2025-12-16
%
% Description:
%   This script demonstrates Active EMI Suppression (AES) on gel phantom data
%   acquired during microwave ablation on a 3T scanner.
%
%
% Inputs (Loaded from .mat):
%   - kdata_rx / kdata_aux  : Raw k-space data (Primary / Auxiliary Coils)
%   - acs_rx / acs_aux      : Auto-calibration signal data (Primary / Auxiliary Coils)
%   - noise_rx / noise_aux  : Noise covariance data (Primary / Auxiliary Coils)
% 
% Image reconstruction results: 
%   - imCC_unc              : Uncorrected 
%   - imCC_decor            : Baseline (noise decorrelation only)
%   - imCC_decor_aes        : Propoased AES (noise decorrelation followed by AES)
% =========================================================================


clear; clc; close all;
restoredefaultpath


%% ---- 1. Path & Environment Setup ----
addpath(genpath('./src'));



%% ---- 2. Load Data & Display Dimensions ----
load('./data/gel_data.mat');
% Get dimensions
[Nkx, Nky, Nch_img] = size(kdata_rx);
[~, ~, Nch_emi] = size(kdata_aux);

% Display dimensions
fprintf('------------------------------------------------\n');
fprintf('Data Loaded Successfully. Dimensions:\n');
fprintf('   > Primary Coil (Nkx, Nky, Ch):   [%d, %d, %d]\n', Nkx, Nky, Nch_img);
fprintf('   > Aux Coil (Nkx, Nky, Ch):       [%d, %d, %d]\n', Nkx, Nky, Nch_emi);
fprintf('------------------------------------------------\n');

%% ---- 3. Processing ----
% A. Noise Decorrelation
fprintf('>> Step 1: Applying Noise Decorrelation... ');
kdata_rx_decor  = apply_decorrelation(kdata_rx, noise_rx);
kdata_aux_decor = apply_decorrelation(kdata_aux, noise_aux);
acs_rx_decor  = apply_decorrelation(acs_rx, noise_rx);
acs_aux_decor  = apply_decorrelation(acs_aux, noise_aux);
fprintf('Done.\n');


% B. Active EMI Suppression (AES)
fprintf('>> Step 2: Performing AES... ');
kdata_rx_decor_aes  = perform_AES(kdata_rx_decor,  kdata_aux_decor);
acs_rx_decor_aes  = perform_AES(acs_rx_decor,  acs_aux_decor);
fprintf('Done.\n');

% C. GRAPPA Reconstruction
fprintf('>> Step 3: Running GRAPPA Reconstruction... ');
% Uncorrected
[im_unc, imCC_unc] = perform_grappa(kdata_rx, acs_rx);
% Baseline
[im_decor, imCC_decor] = perform_grappa(kdata_rx_decor, acs_rx_decor);
% Proposed AES
[im_decor_aes, imCC_decor_aes] = perform_grappa(kdata_rx_decor_aes, acs_rx_decor_aes);
fprintf('Done.\n');
%% ---- 4. Check results----


% Plot 1: k-space Comparison
kdataDisplayRange = [0 1e-4];
shape = [2 6];
kdata_diff = (kdata_rx_decor - kdata_rx_decor_aes);

figure('Name', 'K-Space Comparison', 'Position', [264, 288, 1524, 592]);
tiledlayout(1,3);
nexttile; 
imshow3(abs(kdata_rx_decor(:,1:2:end,:)), kdataDisplayRange, shape);
title('Baseline');
nexttile; 
imshow3(abs(kdata_rx_decor_aes(:,1:2:end,:)), kdataDisplayRange, shape);
title('Proposed AES');
nexttile; 
imshow3(abs(kdata_diff(:,1:2:end,:)), kdataDisplayRange, shape);
title('Complex Difference');


% Plot 2: Single Coil Images
imgDisplayRange = [0 1e-6];
im_diff = (im_decor - im_decor_aes);

figure('Name', 'Single Coil Images', 'Position', [264, 288, 1524, 592]);
tiledlayout(3,1);
nexttile; 
imshow3(imrotate(abs(im_decor), 90), imgDisplayRange, shape);
title('Baseline');
nexttile; 
imshow3(imrotate(abs(im_decor_aes), 90), imgDisplayRange, shape);
title('Proposed AES');
nexttile; 
imshow3(imrotate(abs(im_diff), 90), imgDisplayRange, shape);
title('Complex Difference');


% Plot 3: Coil-Combined Images (Magnitude & Phase)
ccImgDisplayRange = [0 0.7];
shape_cc = [1 1];

figure('Name', 'Coil-Combined Results', 'Position', [264, 288, 1524, 592]);
tiledlayout(2,3);

nexttile;
imshow3(imrotate(abs(imCC_unc), 90), ccImgDisplayRange, shape_cc);
title('Uncorrected - Mag');
nexttile;
imshow3(imrotate(abs(imCC_decor), 90), ccImgDisplayRange, shape_cc);
title('Baseline - Mag');
nexttile;
imshow3(imrotate(abs(imCC_decor_aes), 90), ccImgDisplayRange, shape_cc);
title('Proposed AES - Mag');

nexttile;
imshow3(imrotate(angle(imCC_unc), 90), [], shape_cc);
title('Uncorrected - Phase');
nexttile;
imshow3(imrotate(angle(imCC_decor), 90), [], shape_cc);
title('Baseline - Phase');
nexttile;
imshow3(imrotate(angle(imCC_decor_aes), 90), [], shape_cc);
title('Proposed AES - Phase');