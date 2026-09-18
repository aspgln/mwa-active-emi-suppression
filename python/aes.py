"""
Active EMI Suppression (AES) — Python port of the MATLAB reference.

Implements three steps of the AES pipeline:
  - noise decorrelation (Pruessmann et al., MRM 2001)
  - PCA compression of the auxiliary k-space (SVD over channels)
  - AES itself: build a shifted-aux design matrix and subtract the
    least-squares EMI prediction per primary coil.

The transfer-function modeling for EMI estimation follows the linear
convolution approach in Srinivas et al., MRM 2022 (doi:10.1002/mrm.28992).

This is an idiomatic Python port; results are visually comparable to the
MATLAB pipeline but not bit-identical (different linear-solver behaviour,
no MATLAB-specific Image Processing Toolbox calls). The MATLAB pipeline
in ``src/`` is the reference implementation for the paper.

Please cite:
    Dai Q, et al. Magn Reson Med. 2026. [DOI forthcoming]

Author: Qing Dai (qdai@ucla.edu)
Last modified: 260519
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, solve_triangular


# ----------------------------------------------------------------------
# Noise decorrelation
# ----------------------------------------------------------------------
def apply_decorrelation(kdata: np.ndarray, noise: np.ndarray) -> np.ndarray:
    """Whiten multi-channel k-space using a Cholesky factor of the
    noise covariance estimated from pre-scan noise samples.

    Parameters
    ----------
    kdata : (Nkx, Nky, Nch) complex
        Acquired k-space data.
    noise : (Nkx_pre, Nch, Nky_pre, Navg) complex
        Pre-scan noise samples, same coil ordering as ``kdata``.

    Returns
    -------
    (Nkx, Nky, Nch) complex
        Noise-decorrelated k-space.
    """
    # Flatten noise to (Nch, Nsamples) for covariance estimation.
    nch = noise.shape[1]
    noise_flat = np.moveaxis(noise, 1, 0).reshape(nch, -1)

    # Sample covariance, then normalize so the diagonal has unit mean
    # (matches the MATLAB reference; keeps the Cholesky factor well-scaled).
    psi = noise_flat @ noise_flat.conj().T
    psi = psi / np.mean(np.abs(np.diag(psi)))
    # Force the diagonal to be real & positive before factoring.
    np.fill_diagonal(psi, np.abs(np.diag(psi)))

    # Lower-triangular Cholesky factor of the (Hermitian) covariance.
    L = cholesky(psi, lower=True)

    # Whiten each ky line independently: kdata[:, jj, :].T = L \ x
    nkx, nky, _ = kdata.shape
    out = np.empty_like(kdata)
    for jj in range(nky):
        # L is (Nch, Nch), RHS is (Nch, Nkx) → solve along channel axis.
        rhs = kdata[:, jj, :].T  # (Nch, Nkx)
        out[:, jj, :] = solve_triangular(L, rhs, lower=True).T
    return out


# ----------------------------------------------------------------------
# PCA compression of auxiliary channels
# ----------------------------------------------------------------------
def apply_pca(kdata_aux: np.ndarray, num_components: int = 1):
    """SVD-based coil compression of the auxiliary k-space.

    Returns
    -------
    kdata_pca : (Nkx, Nky, num_components) complex
        Projection of ``kdata_aux`` onto the top principal components.
    energy_spectrum : (Nch_aux,) float
        Normalized variance explained per PC (sums to 1).
    """
    nkx, nky, nch = kdata_aux.shape
    flat = kdata_aux.reshape(nkx * nky, nch)

    # Economy SVD over channels; V holds the right singular vectors as columns.
    _, S, Vh = np.linalg.svd(flat, full_matrices=False)
    V = Vh.conj().T

    energy = S ** 2
    total = energy.sum()
    energy_spectrum = energy / total if total > 0 else np.zeros_like(energy)

    projected = flat @ V[:, :num_components]
    return projected.reshape(nkx, nky, num_components), energy_spectrum


# ----------------------------------------------------------------------
# Acceleration detection (mirrors MATLAB check_accel_ratio along ky)
# ----------------------------------------------------------------------
def check_accel_ratio(kdata: np.ndarray, axis: int = 1):
    """Detect the parallel-imaging acceleration factor and the index of
    the first sampled ky line.

    Looks at a 1-D slice along ``axis`` taken at index 0 of every other
    axis, then counts leading zeros — same heuristic as the MATLAB
    reference.
    """
    # Build the index: first element of every axis except the one we scan.
    idx = [0] * kdata.ndim
    idx[axis] = slice(None)
    patch = kdata[tuple(idx)]

    r_accel = 1
    first_pe = 0  # 0-based here (MATLAB used 1-based)
    if patch[0] == 0:
        # First line is zero → find first non-zero, line 2 is first sampled.
        while patch[r_accel] == 0:
            r_accel += 1
        first_pe = 1
    elif patch.size >= 2 and patch[1] == 0:
        # Second line is zero → find next non-zero starting from index 1.
        while patch[r_accel] == 0:
            r_accel += 1
        first_pe = 0
    return r_accel, first_pe


# ----------------------------------------------------------------------
# AES: predict and subtract EMI from each primary coil
# ----------------------------------------------------------------------
def perform_aes(
    kdata_rx: np.ndarray,
    kdata_aux: np.ndarray,
    num_components: int = 1,
    dkx: int = 2,
    dky: int = 0,
) -> np.ndarray:
    """Active EMI Suppression on multi-channel MRI k-space.

    Parameters
    ----------
    kdata_rx : (Nkx, Nky, Nch_img) complex
        Primary coil k-space.
    kdata_aux : (Nkx, Nky, Nch_emi) complex
        Auxiliary (EMI-sensing) coil k-space.
    num_components : int, optional
        Number of PCA components to retain from the auxiliary channels.
        Default 1. Pass 0 to skip PCA and use all aux channels as-is.
    dkx, dky : int, optional
        Kernel half-widths along readout and phase-encode. Default (2, 0).

    Returns
    -------
    (Nkx, Nky, Nch_img) complex
        EMI-suppressed primary k-space (zero-fills restored on accelerated data).
    """
    nkx, nky, nch_rx = kdata_rx.shape

    r_accel, first_pe = check_accel_ratio(kdata_rx, axis=1)

    if num_components == 0:
        aux = kdata_aux
    else:
        aux, _ = apply_pca(kdata_aux, num_components)

    # Strip zero-filled ky lines for accelerated acquisitions; we'll
    # zero-fill back at the end so the caller sees the original shape.
    if r_accel > 1:
        kdata_rx_used = kdata_rx[:, first_pe::r_accel, :]
        aux_used = aux[:, first_pe::r_accel, :]
    else:
        kdata_rx_used = kdata_rx
        aux_used = aux
    nky_used = kdata_rx_used.shape[1]
    nch_aux = aux_used.shape[2]

    # Build the design matrix E by stacking shifted copies of the aux k-space.
    # Each (col_shift, lin_shift, aux_channel) triple is one column of E.
    aux_padded = np.pad(aux_used, ((dkx, dkx), (dky, dky), (0, 0)))
    n_shifts = (2 * dkx + 1) * (2 * dky + 1) * nch_aux
    E = np.empty((nkx * nky_used, n_shifts), dtype=kdata_rx.dtype)

    col = 0
    for col_shift in range(-dkx, dkx + 1):
        for lin_shift in range(-dky, dky + 1):
            shifted = np.roll(aux_padded, shift=(col_shift, lin_shift), axis=(0, 1))
            cropped = shifted[dkx:dkx + nkx, dky:dky + nky_used, :]  # (Nkx, Nky_used, Nch_aux)
            for j in range(nch_aux):
                E[:, col] = cropped[:, :, j].reshape(-1)
                col += 1

    # Per-coil least-squares solve, then subtract the predicted EMI.
    kdata_clean = np.empty_like(kdata_rx_used)
    for ch in range(nch_rx):
        y = kdata_rx_used[:, :, ch].reshape(-1)
        h, *_ = np.linalg.lstsq(E, y, rcond=None)
        pred = (E @ h).reshape(nkx, nky_used)
        kdata_clean[:, :, ch] = kdata_rx_used[:, :, ch] - pred

    # Restore zero-fills on accelerated data.
    if r_accel > 1:
        out = np.zeros_like(kdata_rx)
        out[:, first_pe::r_accel, :] = kdata_clean
        return out
    return kdata_clean
