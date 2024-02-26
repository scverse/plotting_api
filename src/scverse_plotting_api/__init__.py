from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

from .helpers import plot_decorator

if TYPE_CHECKING:
    from anndata import AnnData
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec


@plot_decorator
def plot_umap(
    fig: Figure, gs: GridSpec, /, adata: AnnData, *, color: str
) -> Figure | None:
    main_ax = fig.add_subplot(gs[0, 0])
    cbar_ax = fig.add_subplot(gs[0, 1])
    x, y = zip(*adata.obsm["X_umap"])
    scatter = main_ax.scatter(
        x, y, s=3, c=sc.get.obs_df(adata, color).values, cmap="viridis"
    )
    fig.colorbar(scatter, cax=cbar_ax)
