function [kdata_pca, energy_spectrum] = apply_pca(kdata_aux, numComponents)
% =========================================================================
%  apply_pca - SVD-based coil compression of the auxiliary k-space channels.
%
%  Inputs:
%    kdata_aux       - Complex k-space data from auxiliary coils
%                      [Nx, Ny, Nch_aux]
%    numComponents   - Scalar integer. Number of principal components to
%                      retain.
%
%  Outputs:
%    kdata_pca       - k-space data projected into the PC basis
%                      [Nx, Ny, numComponents]
%    energy_spectrum - Vector. Normalized energy fraction for each PC.
%
%  Author:  Qing Dai (qdai@ucla.edu)
%  Last modified: 260518
% =========================================================================

    [Nkx, Nky, Nch_aux] = size(kdata_aux);
    n_pixels = Nkx * Nky;

    % svd
    kdata_flat = reshape(kdata_aux, n_pixels, Nch_aux);
    [~, S, V] = svd(kdata_flat, 'econ');

    % Compute Energy Distribution
    singular_values = diag(S);
    energy_raw = singular_values.^2;
    total_energy = sum(energy_raw);

    if total_energy > 0
        energy_spectrum = energy_raw / total_energy;
    else
        energy_spectrum = zeros(size(energy_raw));
    end

    % Project data onto all Principal Components
    kdata_rotated = kdata_flat * V(:, 1:numComponents); 

    % Reshape to original dimensions
    kdata_pca = reshape(kdata_rotated, Nkx, Nky, numComponents);

end