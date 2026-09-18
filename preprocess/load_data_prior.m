% ================================================================
% MRI K-Space Preprocessing Script
%
% Purpose:
%   Convert Siemens .dat raw data files into standardized matrices:
%       - k_img (imaging coils)
%       - k_emi (EMI detection coils)
%       - c_img (imaging coils)
%       - c_emi (EMI detection coils)
%       - noise_img (imaging coils)
%       - noise_emi (EMI detection coils)
%   Each frame is saved separately in one file to support flexible
%   downstream analysis and sharing without raw .dat files.
%
% Author:  Qing Dai
% Created: 2025-12-05
% ================================================================

%% ========== INITIALIZATION =====================================
close all; clear; clc;
restoredefaultpath

% ---- Default Paths ----
SRC_PATH   = '/data2/qingdai/Research/MRMWA_Tmap_v2/';
SAVE_ROOT  = fullfile(SRC_PATH, 'reconData');

% ---- Path Setup ----
addpath(genpath(fullfile(SRC_PATH, 'GRAPPA/')));
addpath(genpath(fullfile(SRC_PATH, 'coil_combine_rapid/')));
addpath(genpath(fullfile(SRC_PATH, 'mapVBVD_files/')));

% ---- Default Experiment Settings ----
DATA_ROOT             = '/data2/qingdai/Data/Animal/';
SCAN_NAME             = '2025_05_08_Pig_Chiang_MWA';   % << change for each dataset
N_EMI_COIL            = 18;                    % number of EMI coils
EC_RANGE              = 3;                     % echo index (can be scalar or array)
T_RANGE_MODE          = 'all';                 % 'all' or e.g. 1:100
CENTER_LINE_N         = 24;                    % # of PE lines for ACS center
DEFAULT_DISPLAY_FRAME = 70;                    % visualization frame index
acs_MODE              = "separated ACS";
 % choose: meas_MID00234_FID55830_2d_tmap_trig_iPAT2_target2_seperateACS_65W_10min.dat

%
% DATA_ROOT             = '/data2/qingdai/Data/Animal/';
% SCAN_NAME             = '2025_01_29_Pig_Lu';   % << change for each dataset
% N_EMI_COIL            = 18;                    % number of EMI coils
% EC_RANGE              = 3;                     % echo index (can be scalar or array)
% T_RANGE_MODE          = 'all';                 % 'all' or e.g. 1:100
% CENTER_LINE_N         = 24;                    % # of PE lines for ACS center
% DEFAULT_DISPLAY_FRAME = 70;                    % visualization frame index

% DATA_ROOT             = '/data2/qingdai/Data/Animal/';
% SCAN_NAME             = '2025_03_13_Pig_Lu';   % << change for each dataset
% N_EMI_COIL            = 18;                    % number of EMI coils
% EC_RANGE              = 3;                     % echo index (can be scalar or array)
% T_RANGE_MODE          = 'all';                 % 'all' or e.g. 1:100
% CENTER_LINE_N         = 24;                    % # of PE lines for ACS center
% DEFAULT_DISPLAY_FRAME = 70;                    % visualization frame index
% % Run target 3... 
%
% 
% DATA_ROOT             = '/data2/qingdai/Data/MR-MWA/';
% SCAN_NAME             = '2024_10_25_MWA_Gel';   % << change for each dataset
% N_EMI_COIL            = 18;                    % number of EMI coils
% EC_RANGE              = 1;                     % echo index (can be scalar or array)
% T_RANGE_MODE          = 'all';                 % 'all' or e.g. 1:100
% CENTER_LINE_N         = 24;                    % # of PE lines for ACS center
% % DEFAULT_DISPLAY_FRAME = 70;                    % visualization frame index
% % choose: meas_MID00109_FID44660_MRT_EXP1_Gel_2min50W.dat


fprintf('>>> Initialized successfully. \n');

%% ========== SELECT & LOAD DATA =================================
cd(fullfile(DATA_ROOT, SCAN_NAME));
[dat_filename, dat_dir] = uigetfile('*.dat', 'Select raw MRI .dat file:');
cd(SRC_PATH);

% ---- Load Raw K-space ----
fprintf('>>> Loading k-space data from %s ...\n', dat_filename);
rparams.decorrelation = 0;
rparams.compressionNum = 0;
dataObj = io.loadCartesianData_v2(dat_dir, dat_filename, rparams);
fprintf('>>> Data loaded successfully.\n');

%% ========== DEFINE SAVE PATH ===================================
% Extract file base name (without .dat)
[~, dat_base, ~] = fileparts(dat_filename);
SAVE_PATH = fullfile(SAVE_ROOT,SCAN_NAME, dat_base);
if ~exist(SAVE_PATH, 'dir')
    mkdir(SAVE_PATH);
end
fprintf('>>> Data will be saved to: %s\n', SAVE_PATH);

%% ========== SEPARATE IMAGE / EMI COILS ==========================
fprintf('>>> Separating imaging and EMI coils ...\n');
N_coils    = size(dataObj.kdata, 3);
N_IMG_COIL = N_coils - N_EMI_COIL;
ec = 3;

[idxImg, idxEmi] = AES_v3.seperate_img_and_aux_coils(...
    dataObj.kdata(:,:,:,ec,end-1), N_IMG_COIL, dataObj.R_accel);

% ---- Quick Visualization (optional) ----
kDisplayRange = [0 1e-4];
figure('Name','Coil Separation Preview');
tiledlayout(1,2);
nexttile; plotter.imshow3(abs(dataObj.kdata(:,dataObj.firstPE:2:end,idxImg,ec,end-2)), kDisplayRange);
title('Imaging Coils');
nexttile; plotter.imshow3(abs(dataObj.kdata(:,dataObj.firstPE:2:end,idxEmi,ec,end-2)), kDisplayRange);
title('EMI Coils');

%% ========== ASSIGN RAW COMPONENTS ==============================
fprintf('>>> Extracting coil-specific data ...\n');
kdata_img = dataObj.kdata(:,:,idxImg,:,:);
kdata_emi = dataObj.kdata(:,:,idxEmi,:,:);
acs_img   = dataObj.acs(:,:,idxImg,:,:);
acs_emi   = dataObj.acs(:,:,idxEmi,:,:);
noise_img  = dataObj.noise(:,idxImg,:,:,:,:);
noise_emi  = dataObj.noise(:,idxEmi,:,:,:,:);

Nt = size(kdata_img, 5); % total number of temporal frames

% ========== HANDLE kdata ECHOES / FRAMES ===============================
tRange = 1:Nt-1;  % remove last frame 

Npe = size(dataObj.kdata,2);
centerLines = util.getCenterLines(Npe, CENTER_LINE_N);

% Select echo and time ranges 
% kdata_img = kdata_img(:,:,:,EC_RANGE,tRange);
% kdata_emi = kdata_emi(:,:,:,EC_RANGE,tRange);

fprintf('>>> Selected %d frame(s) for processing.\n', numel(tRange));

% ========== HANDLE ACS MODES =====================================
% if isempty(acs_img)
%     acs_MODE = "NA";
% elseif size(acs_img,5) > 1
%     acs_MODE = "integrated ACS";
%     acs_img = acs_img(:,centerLines,:,EC_RANGE,tRange);
%     acs_emi = acs_emi(:,centerLines,:,EC_RANGE,tRange);
% else
%     acs_MODE = "separated ACS";
% end
fprintf('>>> ACS mode detected: %s\n', acs_MODE);




%% ========== SAVE FRAME-BY-FRAME (raw k-space data, separated ACS) =================================
fprintf('>>> Saving preprocessed k-space data to:\n%s\n', SAVE_PATH);
SAVE_PATH_kdata = fullfile(SAVE_PATH, 'kdata');
mkdir(SAVE_PATH_kdata);

ec = 3;

for t = 1:numel(tRange)
    k_img = squeeze(kdata_img(:,:,:,ec,t)); % 3D 
    k_emi = squeeze(kdata_emi(:,:,:,ec,t));

    % ---- Save all frame data into ONE file ----
    fname = fullfile(SAVE_PATH_kdata,sprintf('kdata_frame%.3i.mat', tRange(t)));
    save(fname, "k_img", "k_emi" ,  '-v7.3');

    if mod(t,10)==0 || t==numel(tRange)
        fprintf('   Saved frame %d/%d\n', t, numel(tRange));
    end
end


c_img = acs_img;
c_emi = acs_emi;

SAVE_PATH_acsData= fullfile(SAVE_PATH, 'ACS');
mkdir(SAVE_PATH_acsData);
fname = fullfile(SAVE_PATH_acsData,   sprintf('acs_img.mat'));
save(fname, "c_img" , '-v7.3');
fname = fullfile(SAVE_PATH_acsData, sprintf('acs_emi.mat'));
save(fname, "c_emi" ,  '-v7.3');

SAVE_PATH_noiseData= fullfile(SAVE_PATH, 'noise');
mkdir(SAVE_PATH_noiseData);
fname = fullfile(SAVE_PATH_noiseData,   sprintf('noise_img.mat'));
save(fname, "noise_img" , '-v7.3');
fname = fullfile(SAVE_PATH_noiseData, sprintf('noise_emi.mat'));
save(fname, "noise_emi" ,  '-v7.3');
fprintf('   Saved pre-scan noise measurement\n');



fname = fullfile(SAVE_PATH, sprintf('noise.mat'));
save(fname, "noise" , '-v7.3');
fprintf('   Saved pre-scan noise measurement\n');

fprintf('>>> Preprocessing complete! ✅\n');

