from __future__ import annotations

from typing import TYPE_CHECKING, TypedDict, overload

import numpy as np
import scanpy as sc

from scverse_plotting_api.helpers import plot_context

if TYPE_CHECKING:
    from collections.abc import Collection
    from typing import Required, Unpack

    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec, SubplotSpec


class UmapArgs(TypedDict, total=False):
    color: Required[Collection[str]]


@overload
def plot_umap(adata: AnnData, *, ax: None = None, **kw: Unpack[UmapArgs]) -> Figure: ...


@overload
def plot_umap(
    adata: AnnData, *, ax: Axes | GridSpec | SubplotSpec, **kw: Unpack[UmapArgs]
) -> None: ...


def plot_umap(
    adata: AnnData,
    *,
    ax: Axes | GridSpec | SubplotSpec | None = None,
    **kw: Unpack[UmapArgs],
) -> Figure | None:
    with plot_context(ax, n_plots=len(kw["color"])) as (fig, gs):
        _plot_umap_inner(adata, fig=fig, gs=gs, **kw)
    return fig if ax is None else None


def _plot_umap_inner(
    adata: AnnData,
    *,
    fig: Figure,
    gs: GridSpec | GridSpecFromSubplotSpec,
    **kw: Unpack[UmapArgs],
) -> None:
    axs = np.empty(len(kw["color"]), dtype=object)
    x, y = zip(*adata.obsm["X_umap"])
    for i, c in enumerate(kw["color"]):
        axs[i] = fig.add_subplot(gs[i], sharex=axs[0], sharey=axs[0])
        scatter = axs[i].scatter(x, y, s=3, c=sc.get.obs_df(adata, c).values)
    fig.colorbar(scatter, ax=axs)
