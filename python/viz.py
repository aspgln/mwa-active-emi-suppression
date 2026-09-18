"""
Visualization helpers for multi-channel k-space / image data.

Python counterpart of src/utils/imshow3.m (Lustig 2012): tile the channels
of a 3D array into a single montage so all coils can be inspected at once.

Last modified: 260618
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def imshow3(img, vrange=None, shape=None, cmap='gray', title=None, ax=None):
    """Show every channel of a 3D array as a montage (like imshow3.m).

    Parameters
    ----------
    img : (Nx, Ny, Nch) array, real or complex
        Channels live on the last axis. Complex input is shown as magnitude
        (matplotlib can't display complex values directly).
    vrange : (vmin, vmax), optional
        Display window. Defaults to (min, max) of the data — the same default
        as the MATLAB version. Pass e.g. (0, 1e-4) to clip bright k-space DC.
    shape : (nrows, ncols), optional
        Montage grid. Defaults to a near-square layout, padding with blank
        tiles if Nch isn't a perfect fit.
    cmap : str
        Colormap. 'gray' for magnitude, e.g. 'twilight' for phase.
    title : str, optional
        Figure title.
    ax : matplotlib Axes, optional
        Draw into an existing axis instead of making a new figure.

    Returns
    -------
    ax : the axis the montage was drawn on.
    """
    img = np.asarray(img)
    if np.iscomplexobj(img):
        img = np.abs(img)

    nx, ny, nch = img.shape

    # Decide the grid. Near-square by default (mirrors imshow3's ceil(sqrt)).
    if shape is None:
        ncols = int(np.ceil(np.sqrt(nch)))
        nrows = int(np.ceil(nch / ncols))
    else:
        nrows, ncols = shape

    # Pad with blank channels so the grid is exactly filled.
    n_tiles = nrows * ncols
    if n_tiles > nch:
        pad = np.zeros((nx, ny, n_tiles - nch), dtype=img.dtype)
        img = np.concatenate([img, pad], axis=2)

    # Build the montage row by row: hstack tiles in a row, vstack the rows.
    rows = []
    for r in range(nrows):
        tiles = [img[:, :, r * ncols + c] for c in range(ncols)]
        rows.append(np.hstack(tiles))
    montage = np.vstack(rows)

    if vrange is None:
        vrange = (montage.min(), montage.max())

    if ax is None:
        # Size the figure to the montage aspect ratio so tiles stay square-ish.
        fig, ax = plt.subplots(figsize=(2.2 * ncols, 2.2 * nrows))
    ax.imshow(montage, cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    ax.axis('off')
    if title:
        ax.set_title(title)
    return ax
