"""
Plate-reader signal transforms.
"""

import numpy as np


def gaussian_median_filter(
    x: np.ndarray,
    window: int = 5,
    edge: str = "raw",
) -> np.ndarray:
    """
    Modified median filter (Alon et al. 1998, EMBO J. 17:4238).

    For each window of `window` points along the last axis: sort, drop the min and
    max, average the remainder. window=5 averages the central 3 (Alon's canonical
    setting).

    `edge` governs the first/last `window // 2` points, where no full window exists.
    Alon et al. specify no edge rule.
      "raw" (default) — pass the point through unsmoothed. Keeps every timepoint and
        preserves a linear trend exactly. Averaging a clipped window instead would
        shift the endpoint by ~one sample interval, because a clipped window is
        off-centre; on plate-reader ln(OD) that bias runs ~2x the noise it would
        suppress, so smoothing there is worse than not smoothing at all.
      "nan" — emit no value there. No extrapolation, at the cost of those points.

    Output shape == input shape.
    """
    x = np.asarray(x, dtype=float)
    if window < 3:
        return x.copy()
    if edge not in ("raw", "nan"):
        raise ValueError(f"edge must be 'raw' or 'nan'; got {edge!r}")

    half = window // 2
    out = np.full_like(x, np.nan)
    if x.shape[-1] >= window:
        sw = np.lib.stride_tricks.sliding_window_view(x, window, axis=-1)
        out[..., half:-half] = np.sort(sw, axis=-1)[..., 1:-1].mean(axis=-1)
    if edge == "raw":
        out[..., :half] = x[..., :half]
        out[..., -half:] = x[..., -half:]
    return out

def exp_growth_window(
    growth_rate: np.ndarray, threshold: float = 0.95
) -> tuple[int, int]:
    """
    (start_idx, end_idx) of the contiguous region around argmax(growth_rate)
    where growth_rate >= threshold * max(growth_rate). 

    1D input. Walks outward from argmax until the signal first drops below the
    cutoff or hits a NaN.
    """
    growth_rate = np.asarray(growth_rate, dtype=float)
    if not np.isfinite(growth_rate).any():
        raise ValueError("no finite values in growth_rate")
    idx_max = int(np.nanargmax(growth_rate))
    cutoff = threshold * growth_rate[idx_max]
    start = idx_max
    while (
        start > 0
        and np.isfinite(growth_rate[start - 1])
        and growth_rate[start - 1] >= cutoff
    ):
        start -= 1
    end = idx_max
    n = growth_rate.shape[0]
    while (
        end < n - 1
        and np.isfinite(growth_rate[end + 1])
        and growth_rate[end + 1] >= cutoff
    ):
        end += 1
    return start, end


def fit_autofluorescence(
    od: np.ndarray,
    fl: np.ndarray,
    above_floor: np.ndarray,
    order: int = 3,
) -> tuple[np.ndarray, tuple[float, float], np.ndarray]:
    """
    Polynomial of autofluorescence against OD.
    """
    od = np.asarray(od, dtype=float)
    fl = np.asarray(fl, dtype=float)
    keep = (np.isfinite(od) & np.isfinite(fl) & (od > 0)
            & np.asarray(above_floor, dtype=bool))
    if keep.sum() < order + 1:
        raise ValueError(
            f"{int(keep.sum())} usable points for an order-{order} fit; need {order + 1}"
        )
    coeffs = np.polyfit(od[keep], fl[keep], order)
    return coeffs, (float(od[keep].min()), float(od[keep].max())), keep


def autofluorescence_at(coeffs: np.ndarray, od: np.ndarray) -> np.ndarray:
    """Autofluorescence predicted by a curve from fit_autofluorescence."""
    return np.polyval(coeffs, np.asarray(od, dtype=float))
