from __future__ import annotations

import inspect
from contextlib import contextmanager
from functools import wraps
from itertools import islice
from typing import TYPE_CHECKING, TypedDict, overload

from matplotlib import get_backend
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, SubplotSpec

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Generic, ParamSpec, Protocol

    from matplotlib.gridspec import GridSpecFromSubplotSpec

    P = ParamSpec("P")

    class PlottingImpl(Protocol, Generic[P]):
        def __call__(
            self,
            *args: P.args,
            fig: Figure,
            gs: GridSpec | GridSpecFromSubplotSpec,
            **kwds: P.kwargs,
        ) -> None:
            ...

    class PlottingAPI(Protocol, Generic[P]):
        @overload
        def __call__(self, *args: P.args, ax: None = None, **kwds: P.kwargs) -> Figure:
            ...

        @overload
        def __call__(
            self, *args: P.args, ax: Axes | SubplotSpec, **kwds: P.kwargs
        ) -> None:
            ...


def plot_decorator(f: PlottingImpl[P]) -> PlottingAPI[P]:
    @overload
    def wrapper(
        *args: P.args,
        ax: None = None,
        **kwargs: P.kwargs,
    ) -> Figure:
        ...

    @overload
    def wrapper(
        *args: P.args,
        ax: Axes | SubplotSpec,
        **kwargs: P.kwargs,
    ) -> None:
        ...

    @wraps(f)
    def wrapper(
        *args: P.args, ax: Axes | SubplotSpec | None = None, **kwargs: P.kwargs
    ) -> Figure | None:
        with plot_context(ax) as (fig, gs):
            rv = f(*args, fig=fig, gs=gs, **kwargs)  # type: ignore[func-returns-value]
            if rv is not None:
                name = getattr(f, "__name__", "function")
                msg = f"{name}’ definition should not return anything"
                raise TypeError(msg)
        return fig if ax is None else None

    ax_param = inspect.Parameter(
        "ax", inspect.Parameter.POSITIONAL_ONLY, annotation=Axes | SubplotSpec
    )
    del wrapper.__annotations__["fig"]
    del wrapper.__annotations__["gs"]
    wrapper.__annotations__["ax"] = ax_param.annotation

    wrapper.__signature__ = (sig := inspect.signature(f)).replace(  # type: ignore[attr-defined]
        parameters=[
            ax_param,
            *islice(sig.parameters.values(), 2, None),
        ]
    )

    return wrapper


class GSParams(TypedDict):
    nrows: int
    ncols: int
    width_ratios: tuple[float, float]


@contextmanager
def plot_context(
    ax: Axes | SubplotSpec | None = None,
) -> Generator[tuple[Figure, GridSpec | GridSpecFromSubplotSpec], None, None]:
    # prologue
    ax_was_provided = ax is not None
    gridspec_params = GSParams(nrows=1, ncols=2, width_ratios=(0.7, 0.1))
    fig: Figure | None
    gs: GridSpec | GridSpecFromSubplotSpec
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
