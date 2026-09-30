"""Duplicated from gencirc (Conrad Derbyshire: https://github.com/conradderbyshire/gencirc)."""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import FlowCal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


## LOADING ##

@dataclass
class FlowSample:
    # 'None', 'aTc', 'IPTG' or 'calibration'
    condition: str
    file_path: Path
    # Events are read from file_path on demand, not at discovery
    data: Optional[FlowCal.io.FCSData] = None
    calibrated_data: Optional[np.ndarray] = None
    channel_stats: Optional[Dict[str, float]] = None


@dataclass
class ReplicateFlowData:
    calibration: Optional[FlowSample]
    # condition -> samples
    samples: Dict[str, List[FlowSample]]


def load_sample_mapping(filepath: Path) -> dict:
    """
    Sample_ID -> Circuit_ID, from a CSV.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Sample mapping file not found: {filepath}")

    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        raise RuntimeError(
            f"Unexpected error reading sample mapping {filepath}: {e}"
        ) from e

    required_columns = {"Sample_ID", "Circuit_ID"}
    if not required_columns.issubset(df.columns):
        raise ValueError(
            f"Sample mapping file must contain columns {required_columns}, got {df.columns.tolist()}"
        )

    return dict(zip(df["Sample_ID"], df["Circuit_ID"]))


def load_flow_cytometry_data_from_fcs(
    fcs_root_dir: Path,
    calibration_file: Path,
    condition_map: Optional[Dict[str, str]] = None,
) -> ReplicateFlowData:
    """
    Index a replicate's FCS files without reading their event data.

    Samples sit at <fcs_root_dir>/<condition>/<sample>.fcs. condition_map renames the
    directories and skips unmapped ones.
    """
    fcs_root_dir = Path(fcs_root_dir)
    calibration_file = Path(calibration_file)

    if not fcs_root_dir.is_dir():
        raise FileNotFoundError(f"FCS directory does not exist: {fcs_root_dir}")
    if not calibration_file.is_file():
        raise FileNotFoundError(f"Calibration file does not exist: {calibration_file}")

    calibration_sample = FlowSample(
        condition="calibration", file_path=calibration_file
    )
    calibration_resolved = calibration_file.resolve()
    samples: Dict[str, List[FlowSample]] = {}

    for fcs_path in sorted(fcs_root_dir.rglob("*.fcs")):
        if fcs_path.resolve() == calibration_resolved:
            continue

        parts = fcs_path.relative_to(fcs_root_dir).parts
        if len(parts) != 2:
            print(
                f"[WARNING] Skipping {fcs_path}: expected <condition>/<sample>.fcs, got depth {len(parts)}"
            )
            continue

        directory = parts[0]
        if condition_map is None:
            condition = directory
        elif directory in condition_map:
            condition = condition_map[directory]
        else:
            continue

        samples.setdefault(condition, []).append(
            FlowSample(condition=condition, file_path=fcs_path)
        )

    return ReplicateFlowData(calibration=calibration_sample, samples=samples)


## CALIBRATION ##

# Adapted from flowcal_process.py (BCL SharePoint)


def calibration_transforms(
    filename: Path,
    fsc_ssc_list: List[str],
    mefl_values: Optional[List[float]] = None,
    mefl_channels: Optional[List[str]] = None,
    mepe_values: Optional[List[float]] = None,
    mepe_channels: Optional[List[str]] = None,
    gate_fraction: float = 0.3,
    plot: bool = False,
    plot_save_path: Optional[Path] = None,
    plot_format: str = "pdf",
    full_output: bool = False,
):
    """
    MEFL and MEPE calibration transforms from a bead file.

    Returns (to_mefl, to_mepe), plus the fitted standard curve with full_output=True.
    """
    to_mefl = None
    to_mepe = None
    details = None

    b = FlowCal.io.FCSData(str(filename))
    b = FlowCal.transform.to_rfi(b)

    # density2d with full output for gating
    b_g, __, c, _, _ = FlowCal.gate.density2d(
        b, channels=fsc_ssc_list, gate_fraction=gate_fraction, full_output=True
    )

    plot_dir = None
    plot_filename = None
    if plot and plot_save_path is not None:
        plot_dir = str(Path(plot_save_path).parent)
        plot_filename = Path(plot_save_path).stem
        os.makedirs(plot_dir, exist_ok=True)

    if mefl_values is not None and mefl_channels is not None:
        mefl_output = FlowCal.mef.get_transform_fxn(
            b_g,
            mef_values=mefl_values,
            mef_channels=mefl_channels,
            plot=plot,
            plot_dir=plot_dir,
            plot_filename=plot_filename,
            full_output=True,
        )
        to_mefl = mefl_output.transform_fxn
        details = _fit_details(mefl_output)

    if mepe_values is not None and mepe_channels is not None:
        mepe_output = FlowCal.mef.get_transform_fxn(
            b_g,
            mef_values=mepe_values,
            mef_channels=mepe_channels,
            plot=False,
            full_output=True,
        )
        to_mepe = mepe_output.transform_fxn

    if to_mefl is None and to_mepe is None:
        print(
            f"[WARNING] No calibration transform created for {filename}; both MEFL and MEPE inputs were missing."
        )

    if full_output:
        return to_mefl, to_mepe, details
    return to_mefl, to_mepe


def _fit_details(mefl_output) -> dict:
    """Flatten FlowCal's MEFOutput into the fields worth recording per run."""
    params = dict(
        zip(
            mefl_output.fitting["beads_params_names"][0],
            mefl_output.fitting["beads_params"][0],
        )
    )
    all_rfi = np.asarray(mefl_output.statistic["values"][0], dtype=float)
    selected_rfi = np.asarray(mefl_output.selection["rfi"][0], dtype=float)
    selected_mef = np.asarray(mefl_output.selection["mef"][0], dtype=float)

    return {
        "model": mefl_output.fitting["beads_model_str"][0],
        "m": params.get("m"),
        "b": params.get("b"),
        "fl_mef_auto": params.get("fl_mef_auto"),
        "populations_found": len(all_rfi),
        "populations_selected": len(selected_rfi),
        "bead_medians_rfi": ";".join(f"{v:.1f}" for v in all_rfi),
        "bead_medians_rfi_selected": ";".join(f"{v:.1f}" for v in selected_rfi),
        "mefl_values_selected": ";".join(f"{v:.0f}" for v in selected_mef),
    }


def calibrate_cytometry_data(
    fc_data: FlowCal.io.FCSData,
    to_mefl: Optional[Callable[..., Any]] = None,
    mefl_channels: Optional[List[str]] = None,
    to_mepe: Optional[Callable[..., Any]] = None,
    mepe_channels: Optional[List[str]] = None,
) -> FlowCal.io.FCSData:
    """
    Applies MEFL and/or MEPE calibration transforms to FlowCal data.
    """
    calibrated = fc_data
    if to_mefl is not None and mefl_channels is not None:
        calibrated = to_mefl(calibrated, channels=mefl_channels)
    if to_mepe is not None and mepe_channels is not None:
        calibrated = to_mepe(calibrated, channels=mepe_channels)
    return calibrated


## GATING AND PROCESSING ##

def extract_channels(data):
    channels = []
    for c in data.channels:
        if "FSC" not in c and "SSC" not in c and "-H" not in c and "Time" not in c:
            channels.append(c)
    return channels


def load_cytometry_data_au(filename, fsc_ssc_list, gfp_channel, gate_fraction=0.5):
    raw = FlowCal.io.FCSData(str(filename))
    # Only the needed channels to RFI
    sliced_data = FlowCal.transform.to_rfi(
        raw, channels=[fsc_ssc_list[0], fsc_ssc_list[1], gfp_channel]
    )

    # Gate
    s_g = FlowCal.gate.high_low(sliced_data, channels=fsc_ssc_list)
    s_g = FlowCal.gate.density2d(
        s_g,
        channels=fsc_ssc_list,
        bins=1024,
        gate_fraction=gate_fraction,
        xscale="logicle",
        yscale="logicle",
        sigma=10.0,
        full_output=False,
    )
    return s_g, extract_channels(s_g)


def process_single_condition_file(
    file_path: Path,
    fsc_ssc_list: List[str],
    to_mefl: Optional[callable],
    mefl_channels: Optional[List[str]],
    gfp_channel: str,
    drop_nonpositive: bool = True,
) -> Tuple[Optional[np.ndarray], List[str], Dict[str, float]]:
    """
    Gated, calibrated GFP values, the channels, and per-channel medians for one FCS file.
    """
    fc_data, current_file_channels = load_cytometry_data_au(
        file_path, fsc_ssc_list, gfp_channel
    )
    fc_data = calibrate_cytometry_data(
        fc_data, to_mefl=to_mefl, mefl_channels=mefl_channels
    )

    # Extract GFP channel
    gfp_values = None
    if gfp_channel in current_file_channels:
        gfp_values = np.array(fc_data[:, gfp_channel])
        # A log axis cannot draw non-positive values
        if drop_nonpositive:
            gfp_values = gfp_values[gfp_values > 0.0]

    # Collect stats for all channels
    ch_stats = {c: np.median(np.array(fc_data[:, c])) for c in current_file_channels}

    return gfp_values, current_file_channels, ch_stats




## INDEXING ##

def index_replicate(replicate_data, sample_to_circuit):
    """Map (condition, circuit) -> file path, so a circuit's files can be read on demand."""
    index = {}
    for condition, samples in replicate_data.samples.items():
        for sample in samples:
            circuit_id = sample_to_circuit.get(sample.file_path.stem)
            if circuit_id is not None:
                index[(condition, circuit_id)] = sample.file_path
    return index



## PLOTTING ##

# Paul Tol, https://personal.sron.nl/~pault/, as gencirc uses them
CONDITION_COLOURS_FLOW = {
    "None": (119 / 255, 119 / 255, 119 / 255),
    "aTc": (220 / 255, 5 / 255, 12 / 255),
    "IPTG": (25 / 255, 101 / 255, 176 / 255),
    "Both": (0.0, 0.0, 0.0),
}

# None and Both outlined, so the single-inducer fills stay visible
CONDITION_STYLE = {
    "None": dict(fill=False, linestyle="--", linewidth=1.2, alpha=1.0, zorder=-5),
    "aTc": dict(fill=True, linewidth=0, alpha=0.4, zorder=-20),
    "IPTG": dict(fill=True, linewidth=0, alpha=0.4, zorder=-20),
    "Both": dict(fill=False, linewidth=1.2, alpha=1.0, zorder=-10),
}


def density_panel(ax, values_by_condition, colours=CONDITION_COLOURS_FLOW,
                  xlim=(1, 1.5e7), ylim=3):
    """One replicate's MEFL distributions, one line per inducer state."""
    for condition, values in values_by_condition.items():
        sns.kdeplot(x=values, ax=ax, label=condition, color=colours[condition],
                    common_norm=False, cut=0, log_scale=True,
                    **CONDITION_STYLE[condition])
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, ylim)
