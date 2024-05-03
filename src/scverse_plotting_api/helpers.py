from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING, TypedDict

import numpy as np
from matplotlib import get_backend
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, SubplotSpec

if TYPE_CHECKING:
    import sys
    from collections.abc import Generator

    from matplotlib.gridspec import GridSpecFromSubplotSpec


class GSParams(TypedDict, total=False):
    nrows: int
    ncols: int


@contextmanager
def plot_context(
    ax: Axes | SubplotSpec | GridSpec | None = None,
    n_plots: int = 1,
) -> Generator[tuple[Figure, GridSpec | GridSpecFromSubplotSpec], None, None]:
    # prologue
    ax_was_provided = ax is not None

    gridspec_params = infer_gs_params(n_plots)

    fig: Figure | None
    gs: GridSpec | GridSpecFromSubplotSpec
    if not ax_was_provided:
        fig = Figure(layout="constrained")
        gs = GridSpec(figure=fig, **gridspec_params)
    elif isinstance(ax, GridSpec):
        fig = ax.figure
        assert fig is not None
        gs = ax
    elif isinstance(ax, SubplotSpec):
        fig = ax.get_gridspec().figure
        assert fig is not None
        if len(ax.colspan) * len(ax.rowspan) < n_plots:
            gs = ax.subgridspec(**gridspec_params)
        else:
            gs = ax.subgridspec(ncols=len(ax.colspan), nrows=len(ax.rowspan))
    elif isinstance(ax, Axes):
        fig = ax.get_figure()
        assert fig is not None
        if (subplotspec := ax.get_subplotspec()) is None:
            msg = (
                "Ax object must be based on a gridspec, "
                "e.g. generated using plt.subplots()"
            )
            raise ValueError(msg)
        else:
            gs = subplotspec.subgridspec(**gridspec_params)
            ax.remove()
    else:
        msg = (
            f"If specified, ax must be an Axes or SubplotSpec object, not a {type(ax)}"
        )
        raise ValueError(msg)

    yield fig, gs

    # epilogue
    if not ax_was_provided:
        if get_backend() == "module://matplotlib_inline.backend_inline":
            from matplotlib_inline import backend_inline

            backend_inline.show(close=True)
        else:
            fig.show(warn=True)


def infer_gs_params(n_plots: int) -> GSParams:
    from matplotlib import rcParams

    w, h = rcParams["figure.figsize"]
    n_vert, n_horiz = _infer_dimensions(n_plots, w / h)
    return GSParams(nrows=n_vert, ncols=n_horiz)


def _infer_dimensions(n: int, a: float) -> tuple[int, int]:
    widths = np.arange(1, n + 1)
    heights = (n + widths - 1) // widths
    deviations = np.abs((widths / heights) - a)

    index = np.argmin(deviations)
    return heights[index], widths[index]
