from __future__ import annotations

import inspect
from contextlib import contextmanager
from functools import wraps
from itertools import islice
from typing import TYPE_CHECKING, Concatenate, ParamSpec, TypeVar

from matplotlib import get_backend
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, SubplotSpec

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    P = ParamSpec("P")
    R = TypeVar("R")


# TODO: make ax a kwarg. Depends on:
# https://github.com/python/typing/issues/1009
def plot_decorator(
    f: Callable[Concatenate[Figure, GridSpec, P], R],
) -> Callable[Concatenate[Axes | SubplotSpec, P], R]:
    @wraps(f)
    def wrapper(ax: Axes, *args: P.args, **kwargs: P.kwargs) -> R:
        with plot_context(ax) as (fig, gs):
            rv = f(fig, gs, *args, **kwargs)
            if rv is not None:
                msg = f"{f.__name__}’ definition should not return anything"
                raise TypeError(msg)
        if ax is None:
            return fig

    ax_param = inspect.Parameter(
        "ax", inspect.Parameter.POSITIONAL_ONLY, annotation=Axes | SubplotSpec
    )
    del wrapper.__annotations__["fig"]
    del wrapper.__annotations__["gs"]
    wrapper.__annotations__["ax"] = ax_param.annotation

    wrapper.__signature__ = (sig := inspect.signature(f)).replace(
        parameters=[
            ax_param,
            *islice(sig.parameters.values(), 2, None),
        ]
    )

    return wrapper


@contextmanager
def plot_context(
    ax: Axes | SubplotSpec | None = None,
) -> Generator[tuple[Figure, GridSpec], None, None]:
    # prologue
    ax_was_provided = ax is not None
    gridspec_params = {"nrows": 1, "ncols": 2, "width_ratios": (0.7, 0.1)}
    if not ax_was_provided:
        fig = Figure(layout="tight")
        gs = GridSpec(figure=fig, **gridspec_params)
    elif isinstance(ax, SubplotSpec):
        fig = ax.get_gridspec().figure
        assert fig is not None
        gs = ax.subgridspec(**gridspec_params)
    elif isinstance(ax, Axes):
        fig = ax.get_figure()
        assert fig is not None
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

    yield fig, gs

    # epilogue
    if not ax_was_provided:
        if get_backend() == "module://matplotlib_inline.backend_inline":
            from matplotlib_inline import backend_inline

            backend_inline.show(close=True)
        else:
            fig.show(warn=True)
