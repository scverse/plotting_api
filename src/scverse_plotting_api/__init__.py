from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import scanpy as sc

if TYPE_CHECKING:
    from anndata import AnnData


def plot_umap(
    adata: AnnData, color: str, *, ax: plt.Axes | plt.SubplotSpec | None = None
) -> plt.Figure | None:
    # prologue
    ax_was_provided = ax is not None
    gridspec_params = {"nrows": 1, "ncols": 2, "width_ratios": (0.7, 0.1)}
    if not ax_was_provided:
        fig = plt.figure(layout="tight")
        gs = plt.GridSpec(figure=fig, **gridspec_params)
    elif isinstance(ax, plt.SubplotSpec):
        fig = ax.get_gridspec().figure
        gs = ax.subgridspec(**gridspec_params)
    elif isinstance(ax, plt.Axes):
        fig = ax.get_figure()
        if ax.get_subplotspec() is None:
            msg = (
                "Ax object must be based on a gridspec, "
                "e.g. generated using plt.subplots()"
            )
            raise ValueError(msg)
        else:
            gs = ax.get_subplotspec().subgridspec(**gridspec_params)
            ax.remove()
    else:
        msg = "If specified, ax must be an Axes or SubplotSpec object"
        raise ValueError(msg)

    # plotting
    main_ax = fig.add_subplot(gs[0, 0])
    cbar_ax = fig.add_subplot(gs[0, 1])
    x, y = zip(*adata.obsm["X_umap"])
    scatter = main_ax.scatter(
        x, y, s=3, c=sc.get.obs_df(adata, color).values, cmap="viridis"
    )
    fig.colorbar(scatter, cax=cbar_ax)

    # epilogue
    if not ax_was_provided:
        if plt.get_backend() == "module://matplotlib_inline.backend_inline":
            plt.close()
        else:
            plt.show()
        return fig
