function kdata  = apply_decorrelation(kdata, noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This APPLY_DECORRELATION function perform noise decorrelation across the
% multi-channel k-space data. This algorithm was implemented based on
% the research work from Pruessmann et al. MRM 2001 (doi: 10.1002/mrm.1241)
%
% Input:(1) kdata:   This is the acquired multi-channel k-space data with
%                    dimensions [Nkx,Nky,Nch].
%    
%       (2) noise: This is the prescan data for noise decorrelation
%                    across coil channels. It has the dimensions 
%                    [Nkx_prescan,Nch,Nky_prescan,Navg].
%
% Output:(1) kdata: Noise-decorrelated k-space data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdata  = permute(kdata, [1 3 2]);
noise_nws = squeeze(noise());
[col, cha, lin, ave] = size(noise_nws);

% Reshape noise data for covariance calculation (channels x samples)
noisech_nws = zeros(cha, col*lin*ave);
for c = 1:cha
    tmp = squeeze(noise_nws(:,c,:,:));
    noisech_nws(c,:) = tmp(:);
end

% Calculate and normalize noise covariance matrix
psi = noisech_nws*noisech_nws';
R = psi;
R = R./mean(abs(diag(R)));           % Normalize by mean diagonal
R(eye(cha)==1) = abs(diag(R));       % Ensure positive diagonal

% Cholesky decomposition for decorrelation
L = chol(R,'lower');

% Apply decorrelation to k-space data
[Nkx, Nch, Nky] = size(kdata);
kdata = permute(kdata, [2,1,3]); % Put channels first

for jj = 1:Nky
    kdata(:,:,jj) = L\kdata(:,:,jj);
end
kdata = permute(kdata, [2,3,1,4,5]);

end


