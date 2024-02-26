from __future__ import annotations

from typing import TYPE_CHECKING, overload

import scanpy as sc

from scverse_plotting_api.helpers import plot_context

if TYPE_CHECKING:
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.gridspec import SubplotSpec


@overload
def plot_umap(adata: AnnData, color: str, *, ax: None = None) -> Figure:
    ...


@overload
def plot_umap(adata: AnnData, color: str, *, ax: Axes | SubplotSpec) -> None:
    ...


def plot_umap(
    adata: AnnData, color: str, *, ax: Axes | SubplotSpec | None = None
) -> Figure | None:
    with plot_context(ax, ncols=2, width_ratios=(0.7, 0.1)) as (fig, gs):
        main_ax = fig.add_subplot(gs[0, 0])
        cbar_ax = fig.add_subplot(gs[0, 1])
        x, y = zip(*adata.obsm["X_umap"])
        scatter = main_ax.scatter(
            x, y, s=3, c=sc.get.obs_df(adata, color).values, cmap="viridis"
        )
        fig.colorbar(scatter, cax=cbar_ax)
    return fig if ax is None else None
