"""
Reconstruction utilities: image recon (simple FFT or GRAPPA) + coil combine.

These are *not* part of the AES contribution — they're the standard
parallel-imaging steps that happen after AES has cleaned the primary
k-space. We use community implementations (``pygrappa``) so the AES
notebook stays short and idiomatic; the MATLAB reference uses Mark
Chiew's GRAPPA tools and Shu-Fu Shih's adaptive coil-combine toolbox.

Two recon methods are available via ``reconstruct(..., method=...)``:
  - 'simple_fft' : plain centered ifft2, for fully-sampled k-space.
  - 'grappa'     : GRAPPA fills missing lines first, for accelerated data.

Last modified: 260618
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter
from pygrappa import mdgrappa


def _k2img(kdata: np.ndarray) -> np.ndarray:
    """Centered 2D k-space → image transform (matches MATLAB ifftshift/ifft2/fftshift)."""
    return np.fft.fftshift(
        np.fft.ifft2(np.fft.ifftshift(kdata, axes=(0, 1)), axes=(0, 1)),
        axes=(0, 1),
    )


def simple_fft_recon(kdata: np.ndarray) -> np.ndarray:
    """Plain inverse-FFT recon for fully-sampled k-space.

    No parallel-imaging — just transform each coil to image space. Use this
    when k-space has no missing lines (e.g. after taking every other PE line
    so an R=2 set becomes fully sampled).

    Parameters
    ----------
    kdata : (Nkx, Nky, Nch) complex

    Returns
    -------
    img : (Nkx, Nky, Nch) complex
        Coil-resolved image.
    """
    return _k2img(kdata)


def grappa_recon(kdata: np.ndarray, acs: np.ndarray, kernel_size=(5, 5)) -> np.ndarray:
    """Run 2D GRAPPA with auto-padded ACS.

    Parameters
    ----------
    kdata : (Nkx, Nky, Nch) complex
        Undersampled k-space with zero-fills in the skipped lines.
    acs : (Nkx_acs, Nky_acs, Nch) complex
        Auto-calibration signal. Will be zero-padded along readout to
        match ``kdata`` if needed.

    Returns
    -------
    img : (Nkx, Nky, Nch) complex
        Coil-resolved image-space reconstruction (centered ifft2 over kx,ky).

    Notes
    -----
    We use ``pygrappa.mdgrappa`` (multidimensional GRAPPA), NOT ``pygrappa.grappa``.
    In pygrappa 0.26.3 the plain ``grappa`` only fills missing lines on the
    DC-side half of the readout (kx) axis and leaves the far half exactly zero
    (boundary at the DC line). This reproduces even on clean synthetic R=2 data
    and is independent of FFT-shift convention, odd Nky, and PE orientation —
    so it half-fills k-space silently. ``mdgrappa`` fills 100% and works
    directly on our DC-centered data, so no ifftshift/fftshift is needed.
    """
    # Pad ACS along readout to match kdata if shapes differ — matches MATLAB.
    if acs.shape[0] != kdata.shape[0]:
        pad = (kdata.shape[0] - acs.shape[0]) // 2
        acs = np.pad(acs, ((pad, pad), (0, 0), (0, 0)))

    kfull = mdgrappa(kdata, acs, kernel_size=kernel_size, coil_axis=-1)

    return _k2img(kfull)


def rss_coil_combine(img: np.ndarray) -> np.ndarray:
    """Root-sum-of-squares coil combination — MAGNITUDE ONLY.

    Discards phase, so it cannot be used for PRF thermometry. Fast and
    robust; fine for magnitude display / SNR checks.

    Parameters
    ----------
    img : (Nx, Ny, Nch) complex

    Returns
    -------
    (Nx, Ny) float
        Combined magnitude, normalized so max = 1.
    """
    combined = np.sqrt(np.sum(np.abs(img) ** 2, axis=-1))
    return combined / np.max(combined)


def walsh_coil_combine(img: np.ndarray, patch_size=(11, 11)) -> np.ndarray:
    """Walsh adaptive coil combination — PHASE PRESERVING.

    Python port of the MATLAB 'ACC' method (src/coil_combination, ccRapid3D_ACC.m):
    for each pixel, estimate the coil signal covariance over a local patch,
    take its dominant eigenvector as the combination weights w, and combine
    as ``sum(conj(w) * img)``. The combined image is complex, so its phase
    is usable for PRF-shift thermometry (unlike RSS).

    We compute every pixel (no accel/interp shortcut) — fast enough at demo
    sizes via box-filtered covariances + a single batched eigendecomposition.

    Parameters
    ----------
    img : (Nx, Ny, Nch) complex
        Coil-resolved image.
    patch_size : (px, py)
        Local window for covariance estimation. MATLAB default is (11, 11);
        larger = smoother sensitivity maps.

    Returns
    -------
    (Nx, Ny) complex
        Combined image, normalized so max magnitude = 1.
    """
    nx, ny, nch = img.shape

    # Per-pixel coil covariance R[x,y] = local sum of img·imgᴴ over the patch.
    # Box-filter each coil-pair product (scipy's uniform_filter is real-only,
    # so filter real/imag separately). Scaling is irrelevant — we only need
    # the dominant eigenvector direction.
    R = np.empty((nx, ny, nch, nch), dtype=np.complex128)
    for a in range(nch):
        for b in range(nch):
            prod = img[:, :, a] * np.conj(img[:, :, b])
            R[:, :, a, b] = (uniform_filter(prod.real, patch_size, mode='reflect')
                             + 1j * uniform_filter(prod.imag, patch_size, mode='reflect'))

    # Dominant eigenvector per pixel. R is Hermitian PSD, so eigh gives
    # ascending real eigenvalues; the last column is the dominant direction.
    _, evecs = np.linalg.eigh(R)
    w = evecs[:, :, :, -1]                       # (Nx, Ny, Nch)

    # eigh leaves an arbitrary per-pixel phase on w, which would corrupt the
    # combined phase map. Reference it to the highest-SNR coil (same coil for
    # every pixel) so the phase is smooth and deterministic.
    ref = int(np.argmax(np.sum(np.abs(img) ** 2, axis=(0, 1))))
    w = w * np.exp(-1j * np.angle(w[:, :, ref:ref + 1]))

    combined = np.sum(np.conj(w) * img, axis=-1)
    return combined / np.max(np.abs(combined))


def reconstruct(kdata: np.ndarray, acs: np.ndarray = None,
                method: str = 'simple_fft', combine: str = 'adaptive'):
    """Image recon + coil combine, with selectable recon and combine methods.

    Parameters
    ----------
    kdata : (Nkx, Nky, Nch) complex
        Primary-coil k-space.
    acs : (Nkx_acs, Nky_acs, Nch) complex, optional
        Auto-calibration signal. Required only for method='grappa'.
    method : {'simple_fft', 'grappa'}
        'simple_fft' — plain ifft2 (fully-sampled data); ignores ``acs``.
        'grappa'     — GRAPPA fills missing lines using ``acs`` first.
    combine : {'adaptive', 'rss'}
        'adaptive' — Walsh adaptive combine, complex output (phase preserved,
                     needed for PRF thermometry). Default.
        'rss'      — root-sum-of-squares, magnitude only (no phase).

    Returns
    -------
    img : (Nx, Ny, Nch) complex
        Coil-resolved image.
    img_cc : (Nx, Ny)
        Coil-combined image (complex if 'adaptive', real if 'rss').
    """
    if method == 'simple_fft':
        img = simple_fft_recon(kdata)
    elif method == 'grappa':
        if acs is None:
            raise ValueError("method='grappa' needs an ACS (acs=...).")
        img = grappa_recon(kdata, acs)
    else:
        raise ValueError(f"unknown method '{method}' (use 'simple_fft' or 'grappa')")

    if combine == 'adaptive':
        img_cc = walsh_coil_combine(img)
    elif combine == 'rss':
        img_cc = rss_coil_combine(img)
    else:
        raise ValueError(f"unknown combine '{combine}' (use 'adaptive' or 'rss')")

    return img, img_cc
