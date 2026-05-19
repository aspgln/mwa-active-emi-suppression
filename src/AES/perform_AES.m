function kdata_clean = perform_AES(kdata_rx, kdata_aux, numComponents, dkx, dky)
% =========================================================================
%  perform_AES - Active EMI Suppression (AES) for MRI-Guided Microwave
%  Ablation.
%
%  The transfer function modeling for EMI estimation follows the linear
%  convolution approach described in Srinivas et al., MRM 2022
%  (doi:10.1002/mrm.28992).
%
%  Please cite:
%    Dai Q, et al. Magn Reson Med. 2026. [DOI forthcoming]
%
%  Inputs:
%    kdata_rx      - Primary coil k-space data    [Nkx, Nky, Nch_img]
%    kdata_aux     - Auxiliary coil k-space data  [Nkx, Nky, Nch_emi]
%    numComponents - (optional) Number of PCA components to retain from the
%                    auxiliary channels. Default = 1.
%                    Use 0 to skip PCA and pass aux channels through as-is.
%    dkx           - (optional) Kernel half-width along readout (kx).
%                    Default = 2.
%    dky           - (optional) Kernel half-width along phase-encode (ky).
%                    Default = 0.
%
%  Output:
%    kdata_clean   - EMI-suppressed primary coil k-space data
%
%  Author:  Qing Dai (qdai@ucla.edu)
%  Created: 2025-12-16
%  Last modified: 260518
%  License: See LICENSE
% =========================================================================

% --- Defaults for optional arguments ---
if nargin < 3
    numComponents = 1;
end

if nargin < 5
    dkx = 2;
    dky = 0;
end

% -------------------------------------------------------------------------
% Preparation and Pre-processing
% -------------------------------------------------------------------------
[Nkx, Nky, Nch_rx] = size(kdata_rx);

% Detect acceleration factor (Parallel Imaging)
[R_accel, firstPE] = check_accel_ratio(kdata_rx, 2);


if numComponents == 0
  % Skip PCA entirely — pass all aux channels directly
  kdata_aux_in = kdata_aux;
else
  % Apply PCA compression to aux channels
  [kdata_aux_in, ~] = apply_pca(kdata_aux, numComponents);
end


% Remove zero-filled lines if data is accelerated (undersampled)
if R_accel > 1
    kdata_rx = kdata_rx(:, firstPE:R_accel:end, :);
    kdata_aux_in = kdata_aux_in(:, firstPE:R_accel:end, :);
    Nky2 = size(kdata_rx, 2);
    disp('skipped lines....')
else
    Nky2 = Nky;
end 

% -------------------------------------------------------------------------
% Convolution Matrix Construction
% -------------------------------------------------------------------------
df = zeros(Nkx + 2*dkx, Nky2 + 2*dky, size(kdata_aux_in,3));
for j = 1:size(kdata_aux_in,3)
     df(:,:,j) = padarray(kdata_aux_in(:,:,j), [dkx, dky]);
end

% Build the design matrix 'conv_mat'.
conv_mat = [];
for col_shift = -dkx:dkx
    for lin_shift = -dky:dky
        for j = 1:size(kdata_aux_in,3)
            dftmp = circshift(df(:,:,j), [col_shift, lin_shift]);
            conv_mat = cat(3, conv_mat, ...
                dftmp(dkx+1:end-dkx, dky+1:end-dky));
        end
    end
end

% Reshape into 2D matrix.
E = reshape(conv_mat, size(conv_mat,1)*size(conv_mat,2), size(conv_mat,3));

% -------------------------------------------------------------------------
% EMI Prediction and Suppression
% -------------------------------------------------------------------------
kdata_clean = zeros(Nkx, Nky2, Nch_rx);

for ch = 1:Nch_rx
    kdata_rx_ch = kdata_rx(:, :, ch);
    
    % Solve linear system
    h(:,ch) = E \ kdata_rx_ch(:);
    
    % Synthesize the estimated EMI
    pred = E * h(:,ch);
    pred = reshape(pred, Nkx, Nky2); 
    
    % Subtract estimated EMI from original data
    kdata_clean(:,:,ch) = kdata_rx_ch - pred; 
end 

% Restore zero-filling
if R_accel > 1
    tmp = kdata_clean;
    kdata_clean = zeros(Nkx, Nky, Nch_rx);
    kdata_clean(:, firstPE:R_accel:end, :) = tmp;
end 

end