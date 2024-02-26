from __future__ import annotations

from typing import TYPE_CHECKING, overload

import scanpy as sc

from .helpers import plot_context, plot_decorator

if TYPE_CHECKING:
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec


@overload
def plot_umap(ax: Axes, adata: AnnData, *, color: str) -> Figure:
    ...


@overload
def plot_umap(ax: None, adata: AnnData, *, color: str) -> None:
    ...


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


@overload
def plot_umap_manual(adata: AnnData, *, color: str, ax: None = None) -> Figure:
    ...


@overload
def plot_umap_manual(adata: AnnData, *, color: str, ax: Axes) -> None:
    ...


def plot_umap_manual(
    adata: AnnData, *, color: str, ax: Axes | None = None
) -> Figure | None:
    with plot_context(ax) as (fig, gs):
        main_ax = fig.add_subplot(gs[0, 0])
        cbar_ax = fig.add_subplot(gs[0, 1])
        x, y = zip(*adata.obsm["X_umap"])
        scatter = main_ax.scatter(
            x, y, s=3, c=sc.get.obs_df(adata, color).values, cmap="viridis"
        )
        fig.colorbar(scatter, cax=cbar_ax)
    if ax is None:
        return fig
