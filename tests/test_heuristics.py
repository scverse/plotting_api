from __future__ import annotations

import pytest

from scverse_plotting_api.helpers import _infer_dimensions


@pytest.mark.parametrize(
    ("figsize", "n_plots", "rowbycol"),
    [
        pytest.param((12.0, 3.0), 3, (1, 3), id="wide"),
        pytest.param((3.5, 12.0), 3, (3, 1), id="tall"),
        pytest.param((3.5, 4.0), 8, (3, 3), id="square"),
    ],
)
def test_infer_dimensions(
    figsize: tuple[float, float], n_plots: int, rowbycol: tuple[int, int]
) -> None:
    assert _infer_dimensions(n_plots, figsize[0] / figsize[1]) == rowbycol
