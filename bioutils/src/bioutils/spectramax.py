"""SpectraMax plate reader file parsing.

File format: UTF-16 encoded TSV with multiple tables at different line ranges.
"""

import logging
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)


def read_plate_file(filepath: str) -> list[str]:
    """Read a SpectraMax plate reader file (UTF-16 encoded)."""
    with open(filepath, "r", encoding="utf-16") as f:
        return f.readlines()


def reshape(data: np.ndarray, nrows: int, ncols: int) -> np.ndarray:
    """Reshape flat well data [time, wells] to plate format [row, col, time]."""
    nsamples, nwells = data.shape
    assert nwells == nrows * ncols, f"expected {nrows * ncols} wells, got {nwells}"
    return data.T.reshape(nrows, ncols, nsamples)


def _parse_timestamp_seconds(timestamp: str) -> float:
    t = datetime.strptime(timestamp, "%H:%M:%S").time()
    return t.hour * 3600 + t.minute * 60 + t.second


def parse_tsv_table(lines: list[str], n_wells: int = 384) -> dict:
    """
    Parse a single TSV table from plate reader output.

    Expected format (first 2 lines are metadata/headers, skipped):
        Time    Temperature(°C)    A1    A2    ...
        00:00:00    37    0.0946    0.092    ...
    """
    time_seconds: list[float] = []
    temp: list[float] = []
    well_data: list[list[float]] = []

    for i, line in enumerate(lines[2:]):
        if not line.strip():
            continue
        values = line.strip().split("\t")
        time_seconds.append(_parse_timestamp_seconds(values[0]))
        temp.append(float(values[1]))
        try:
            well_data.append([float(v) for v in values[2:]])
        except ValueError as e:
            logger.warning("line %d: failed to parse well values: %s", i, e)
            continue

    logger.info("loaded %d timepoints", len(time_seconds))
    return {
        "time": np.array(time_seconds),
        "temp": np.array(temp),
        "data": np.array(well_data),
    }


def time_to_minutes(time_seconds: np.ndarray) -> np.ndarray:
    """Convert seconds to minutes."""
    return time_seconds / 60.0
