"""Configuration shared by the analyse_pr notebooks.
"""

from pathlib import Path

import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter

from bioutils import spectramax
from bioutils.microplate import get_adjacent_wells, get_circuit_wells
from bioutils.transform import gaussian_median_filter


## PATHS ##
REPO_DIR = Path(__file__).resolve().parent.parent

# TODO - To run for replicates CD1, CD2 download the data from https://zenodo.org/records/18378852
SUPPLEMENTARY_DIR = REPO_DIR / "data/masters_thesis_supplementary"
GENCIRC_DIR = SUPPLEMENTARY_DIR / "data/gencirc"


## SMOOTHING ##
# (Alon et al. 1998, doi:10.1093/emboj/17.15.4238)
# hybrid Gaussian median filter; window 5 averages the central 3
GMF_WINDOW = 5
def smooth(series):
    """
    GMF Smoothing on all raw data
    (Alon et al. 1998, doi:10.1093/emboj/17.15.4238)
    """
    return gaussian_median_filter(
        np.asarray(series, dtype=float),
        window=GMF_WINDOW, edge="raw")


## FIGURES ##

# 100 for WIP, 900 for final
FIGURE_DPI = 100
PDF_DPI = 72
# keep equal to seq/sbol.py DPI
PNG_DPI = 300

# Inducer conditions in cond_labels order, so the calibration panel, the method scatter and the
# per-circuit pages colour the same condition alike.
CONDITION_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red")

REPLICATE_MARKERS = ("x", "o", "^")

LABEL_FONTSIZE = 9
TICK_FONTSIZE = 8
LEGEND_FONTSIZE = 8


## PLATES ##

# Every circuit's wells; identical on all three plates.
CIRCUIT_WELLS = {circuit: get_circuit_wells(circuit) for circuit in range(1, 93)}

# The three screening plates.
REPLICATES = [
    {
        "tag": "Rep1",
        "label": "CD1",
        "date": "2025-05-13",
        "file": GENCIRC_DIR / "raw/Michiel Library/25_05_13_ML_384_well_experiment.txt",
        # 1-indexed, inclusive
        "od_table_lines": (2, 44),
        "fl_table_lines": (47, 89),
        "wells": {
            "circuits": CIRCUIT_WELLS,
            "pRPU": get_adjacent_wells(14, 16, 4),
            "pL2": get_adjacent_wells(15, 16, 4),
            "wt": get_adjacent_wells(14, 20, 4),
            "blank": get_adjacent_wells(15, 20, 4),
        },
        "cond_labels": ("None", "aTc", "IPTG", "Both"),
    },
    {
        "tag": "Rep2",
        "label": "CD2",
        "date": "2025-05-28",
        "file": GENCIRC_DIR / "raw/Michiel Library/25_05_28_ML_384_well_experiment.txt",
        "od_table_lines": (2, 44),
        "fl_table_lines": (47, 89),
        "wells": {
            "circuits": CIRCUIT_WELLS,
            "pRPU": get_adjacent_wells(14, 16, 4),
            "pL2": get_adjacent_wells(15, 16, 4),
            "wt": get_adjacent_wells(14, 20, 4),
            "blank": get_adjacent_wells(15, 20, 4),
        },
        "cond_labels": ("None", "aTc", "IPTG", "Both"),
    },
    {
        "tag": "Rep3",
        "label": "AL1",
        "date": "2025-10-21",
        # dataset included in this repo
        "file": REPO_DIR / "data/plate_reader/2025_10_21_ML3_384/conrad_25_05_05_ML_384_well_experiment.txt",
        "od_table_lines": (2, 44),
        "fl_table_lines": (47, 89),
        "wells": {
            "circuits": CIRCUIT_WELLS,
            "pRPU": get_adjacent_wells(14, 16, 4),
            "pL2": get_adjacent_wells(15, 16, 4),
            "wt": get_adjacent_wells(14, 20, 4),
            "blank": get_adjacent_wells(15, 20, 4),
        },
        "cond_labels": ("None", "aTc", "IPTG", "Both"),
    },
]

TAGS = [rep["tag"] for rep in REPLICATES]

# Plate layout is identical across replicates, so wells resolve once.
WELLS = REPLICATES[0]["wells"]


def load_plate_series():
    """Raw OD, raw fluorescence and the time grid for every replicate, keyed by tag.
    """
    od_raw, fl_raw, t_minutes = {}, {}, {}
    for rep in REPLICATES:
        lines = spectramax.read_plate_file(rep["file"])
        od_start, od_end = rep["od_table_lines"]
        fl_start, fl_end = rep["fl_table_lines"]
        od_table = spectramax.parse_tsv_table(lines[od_start - 1:od_end])
        fl_table = spectramax.parse_tsv_table(lines[fl_start - 1:fl_end])
        od_raw[rep["tag"]] = spectramax.reshape(od_table["data"], 16, 24)
        fl_raw[rep["tag"]] = spectramax.reshape(fl_table["data"], 16, 24)
        t_minutes[rep["tag"]] = spectramax.time_to_minutes(od_table["time"])

    reference = t_minutes[TAGS[0]]
    for tag in TAGS[1:]:
        if (t_minutes[tag].shape != reference.shape
                or not np.allclose(t_minutes[tag], reference, atol=1.0)):
            raise ValueError(
                f"t_minutes differs between {TAGS[0]} and {tag} by more than 1 min")
    return od_raw, fl_raw, t_minutes


## LOW SATURATION FLOOR ##

# (Beal et al. 2020, doi:10.1038/s42003-020-01127-5)
# "low saturation by indistinguishability from the blank replicates, here set to be anything
# <2 blank standard deviations above the mean blank"
BLANK_SD_MULTIPLE = 2



## CONTAMINATED BLANK ##

# AL1's None-condition blank well (P21) grew during the run.
# It is dropped from the pool.
CONTAMINATED_BLANK_PLATE = "AL1"
BLANK_EXCLUDE = {CONTAMINATED_BLANK_PLATE: 0}


def usable_blank_wells(n_wells, plate):
    """Blank well indices for a plate, with the contaminated one dropped.
    """
    return [well for well in range(n_wells) if well != BLANK_EXCLUDE.get(plate, -1)]


def pooled_blank_and_threshold(blank_wells, plate):
    """One plate's smoothed blank series and per-timepoint detection threshold.
    """
    usable = usable_blank_wells(np.shape(blank_wells)[0], plate)
    smoothed = smooth(blank_wells)[usable]
    return (smoothed.mean(axis=0),
            BLANK_SD_MULTIPLE * smoothed.std(axis=0, ddof=1))


def pooled_blank_fl(blank_wells, plate):
    """One plate's blank fluorescence, pooled over its usable wells.
    """
    return smooth(blank_wells)[usable_blank_wells(np.shape(blank_wells)[0], plate)].mean(axis=0)


def biomass_and_floor(raw_wells, pooled_blank, threshold):
    """Smoothed blank-corrected biomass, and the detection mask tested on it.
    """
    corrected = smooth(raw_wells) - pooled_blank
    return corrected, corrected >= threshold


def blank_series(od_raw, fl_raw):
    """Per-plate pooled blank, detection threshold and pooled blank fluorescence.
    """
    pooled_blank, threshold, pooled_fluorescence = {}, {}, {}
    for rep in REPLICATES:
        tag, plate = rep["tag"], rep["label"]
        rows, cols = zip(*rep["wells"]["blank"])
        pooled_blank[tag], threshold[tag] = pooled_blank_and_threshold(
            od_raw[tag][list(rows), list(cols), :], plate)
        pooled_fluorescence[tag] = pooled_blank_fl(
            fl_raw[tag][list(rows), list(cols), :], plate)
    return pooled_blank, threshold, pooled_fluorescence


## SAMPLING POINT ##
SAMPLING_METHOD = "sg"


## AUTOFLUORESCENCE CALIBRATION ##

# Timepoints [0:cutoff] of each pL2 well feed the calibration fit.
# Chosen by visual analysis.
PL2_FIT_CUTOFF = {
    "CD1": {"None": 22, "aTc": 22, "IPTG": 22, "Both": 22},
    "CD2": {"None": 22, "aTc": 22, "IPTG": 22, "Both": 22},
    "AL1": {"None": 22, "aTc": 22, "IPTG": 22, "Both": 22},
}

# Polynomial order for calibration curve fit (F vs OD).
AUTOFLUORESCENCE_ORDER = 3



## FLOW CYTOMETRY ##
# seq/flow_invalid.txt contains circuits that showed clear multimodal distributions
FLOW_INVALID_FILE = REPO_DIR / "seq/flow_invalid.txt"
FLOW_INVALID = {f"Circuit {int(line)}" for line in FLOW_INVALID_FILE.read_text().split()}


## GROWTH RATE ##

# Fraction of the max growth rate taken as the exponential phase
# (GrowthRates, doi:10.1093/molbev/mst187)
EXP_GROWTH_PHASE_THRESHOLD = 0.95

# Polyfit
# (Ronen et al. 2002, doi:10.1073/pnas.152046799)
# source fits the sixth-order polynomial to log promoter activity; here to ln(OD)
POLY_ORDER = 6
# Excluded edge points
POLY_EDGE_EXCLUDE = 4


# Savitzky-Golay
# (Vidal et al. 2022, doi:10.1093/synbio/ysac020)
SG_POLYORDER = 3
SG_WINDOW = 11
# Excluded edge points
SG_EDGE_EXCLUDE = SG_WINDOW // 2



## SYNTHESIS RATE ##

# Savitzky-Golay on G, for dG/dt.
G_SG_WINDOW = 7
G_SG_POLYORDER = 3
G_SG_EDGE_EXCLUDE = (G_SG_WINDOW - 1) // 2


def dg_dt(g_series, t_minutes):
    """Savitzky-Golay first derivative of G.
    """
    g = np.atleast_2d(np.asarray(g_series, dtype=float))
    t = np.asarray(t_minutes, dtype=float)
    half = G_SG_EDGE_EXCLUDE
    derivative = np.full(g.shape, np.nan)
    for index in range(half, g.shape[-1] - half):
        window = g[:, index - half:index + half + 1]
        complete = np.isfinite(window).all(axis=-1)
        if not complete.any():
            continue
        offsets = t[index - half:index + half + 1] - t[index]
        design = np.vander(offsets, G_SG_POLYORDER + 1, increasing=True)
        # row 1 of the pseudo-inverse is the linear coefficient: the slope at offset 0
        derivative[complete, index] = np.linalg.pinv(design)[1] @ window[complete].T
    return derivative.reshape(np.shape(g_series))


def replicate_mean(matrix, wt_floor, pl2_band):
    """Arithmetic mean down the replicate axis.
    """
    kept = np.where(matrix >= wt_floor, matrix, np.nan)
    counted = np.isfinite(kept).sum(axis=0)
    raw = np.where(counted > 0,
                   np.nansum(kept, axis=0) / np.where(counted > 0, counted, 1), np.nan)

    # push values below zero up to the L2 band
    imputed = np.isfinite(raw) & (raw <= 0)
    return np.where(imputed, (pl2_band[0] + pl2_band[1]) / 2, raw), imputed, raw


## RPU PLOT ##
# RPU is plotted on linear-log axis, with squished linear region for small values

# Keep the linear region of RPU plots small
RPU_LINSCALE = 0.2

# Clip very negative values
BARS_PNG_YMIN = -0.5


def rpu_bands(synthesis_rates, rpu_column):
    """(WT band, pL2 band, their union, the axis limits) in RPU, from the synthesis-rate table.
    """
    def zero_band(group):
        values = synthesis_rates.loc[synthesis_rates["Group"] == group, rpu_column]
        return min(0.0, float(values.min())), max(0.0, float(values.max()))

    wt_band, pl2_band = zero_band("WT"), zero_band("pL2")
    band = (min(wt_band[0], pl2_band[0]), max(wt_band[1], pl2_band[1]))
    # pRPU is the reference and reads 1.0 by construction, so it sets no range
    finite = synthesis_rates.loc[synthesis_rates["Group"] != "pRPU", rpu_column].dropna()
    return wt_band, pl2_band, band, (min(finite.min() * 1.5, band[0]), finite.max() * 1.5)


def set_rpu_scale(ax, band, limits):
    """Apply the band-linear y scale and its decade ticks.
    """
    low, high = band

    def forward(value):
        """RPU -> axis position: linear inside the band, log outside it.
        """
        value = np.asarray(value, dtype=float)
        position = np.empty_like(value)
        above, below = value > high, value < low
        inside = ~(above | below)
        position[inside] = RPU_LINSCALE * (value[inside] - low) / (high - low)
        with np.errstate(divide="ignore", invalid="ignore"):
            position[above] = RPU_LINSCALE + np.log10(value[above] / high)
            position[below] = -np.log10(value[below] / low)
        return position

    def inverse(position):
        """Axis position -> RPU.
        """
        position = np.asarray(position, dtype=float)
        value = np.empty_like(position)
        above, below = position > RPU_LINSCALE, position < 0
        inside = ~(above | below)
        value[inside] = low + (high - low) * position[inside] / RPU_LINSCALE
        value[above] = high * 10.0 ** (position[above] - RPU_LINSCALE)
        value[below] = low * 10.0 ** (-position[below])
        return value

    major, minor = [0.0], []
    for start, stop, sign in ((high, limits[1], 1.0), (-low, -limits[0], -1.0)):
        for power in range(int(np.floor(np.log10(start))), int(np.ceil(np.log10(stop))) + 1):
            step = 10.0 ** power
            if start <= step <= stop:
                major.append(sign * step)
            minor.extend(sign * step * k for k in range(2, 10) if start <= step * k <= stop)
    ax.set_yscale("function", functions=(forward, inverse))
    ax.yaxis.set_major_locator(FixedLocator(sorted(major)))
    ax.yaxis.set_minor_locator(FixedLocator(sorted(minor)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
