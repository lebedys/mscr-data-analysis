"""SBOL diagrams of genotypes."""

import io
from pathlib import Path

import matplotlib.pyplot as plt
import parasbolv
from matplotlib.patches import Polygon
from PIL import Image

HERE = Path(__file__).resolve().parent

# glyphs for sbol diagrams
GLYPH_DIR = HERE / "glyphs"

# appended to every construct
DOWNSTREAM = ("GFP", "tEND")

GLYPH = {
    "AraC": "CDS",
    "PhlF": "CDS",
    "SrpR": "CDS",
    "BM3R1": "CDS",
    "TraRW": "CDS",
    "GFP": "CDS",
    "pTac": "Promoter",
    "pTet": "Promoter",
    "pSrpR": "Promoter",
    "pBAD": "Promoter",
    "pTra": "Promoter",
    "pBM3R1": "Promoter",
    "tSrpR": "Terminator",
    "tPhlF": "Terminator",
    "tEND": "Terminator",
    "Spacer1": "InertDNASpacer",
    "Spacer2": "InertDNASpacer",
    "Spacer3": "InertDNASpacer",
}

# a repressor's promoter and terminator carry its colour; the two inducible inputs share grey
COLOUR = {
    "AraC": "#f1932d",
    "PhlF": "#b178a6",
    "SrpR": "#7bb0de",
    "BM3R1": "#5289c7",
    "TraRW": "#d6c1de",
    "GFP": "#4ebb65",
    "pBAD": "#f1932d",
    "pSrpR": "#7bb0de",
    "pBM3R1": "#5289c7",
    "pTra": "#d6c1de",
    "pTac": "#777777",
    "pTet": "#777777",
    "tSrpR": "#7bb0de",
    "tPhlF": "#b178a6",
}

# the SVG path ids a part's colour reaches, and whether it fills or outlines them
STYLE_PATHS = {
    "CDS": [("cds", "facecolor")],
    "Promoter": [("promoter-body", "edgecolor"), ("promoter-head", "edgecolor")],
    "Terminator": [("terminator-body", "edgecolor"), ("terminator-head", "edgecolor")],
    "Operator": [("operator-path", "edgecolor")],
}

# inducible terminators
INDUCIBLE = {"tSrpR", "tPhlF"}

RULES = [
    ("AraC", "pBAD", "stimulation"),
    ("SrpR", "pSrpR", "inhibition"),
    ("BM3R1", "pBM3R1", "inhibition"),
    ("SrpR", "tSrpR", "stimulation"),
    ("PhlF", "tPhlF", "stimulation"),
]

# offset from each glyph's centroid, tuned so every label outside a CDS lands in one row
LABEL_SKEW = {
    "CDS": (0, 0),
    "Promoter": (0, -18),
    "Terminator": (0, -15),
    "InertDNASpacer": (0, -8),
}

# the glyph already says what these are
UNLABELLED = {"Spacer1", "Spacer2", "Spacer3", "tEND"}

LABEL_PT = 5

# spacing between glyphs, in the glyphs' own units
GAP = 3.0

# whitespace around the construct, and above the topmost arc
MARGIN = 4.0
MARGIN_ABOVE = 6.0

# gap between a glyph and the arc that leaves or lands on it
CLEARANCE = 5.0
HEAD_WIDTH = 5.0
HEAD_HEIGHT = 5.0

# shortest vertical that still reads as a descender rather than as a head sitting on the line
DESCENT = 8.0
LEVEL_STEP = 9.0

# sizing
CLASSIFY_WIDTH_IN = 8.5
BARS_WIDTH_IN = 5.5

# a whole multiple of FIGURE_DPI
DPI = 300

# parsed once on import rather than once per circuit
RENDERER = parasbolv.GlyphRenderer(glyph_path=str(GLYPH_DIR))


def draw(ax, parts, width_in):
    """Glyphs left to right over a baseline; returns each drawn site and the content bounds."""
    x = 0.0
    xs, ys, label_ys, sites = [], [], [], []
    for part in parts:
        glyph = GLYPH[part]
        operator = None
        # the operator takes the terminator's own colour, so the pair reads as one site
        if part in INDUCIBLE:
            operator, (x, _) = RENDERER.draw_glyph(
                ax, "Operator", (x, 0),
                user_parameters={"label_parameters": None},
                user_style={path: {key: COLOUR[part]}
                            for path, key in STYLE_PATHS["Operator"]},
            )
            xs += [operator[0][0], operator[1][0]]
            ys += [operator[0][1], operator[1][1]]
            x += GAP
        label = None if part in UNLABELLED else {
            "text": part,
            "color": (0, 0, 0),
            "xy_skew": LABEL_SKEW[glyph],
            "userfont": {"size": LABEL_PT},
        }
        style = None if part not in COLOUR else {
            path: {key: COLOUR[part]} for path, key in STYLE_PATHS[glyph]
        }
        bounds, (x, _) = RENDERER.draw_glyph(
            ax, glyph, (x, 0),
            user_parameters={"label_parameters": label},
            user_style=style,
        )
        xs += [bounds[0][0], bounds[1][0]]
        ys += [bounds[0][1], bounds[1][1]]
        if label is not None:
            label_ys.append(ax.texts[-1].get_position()[1])
        # an arc targets the operator and its terminator together
        if operator is not None:
            bounds = ((operator[0][0], min(operator[0][1], bounds[0][1])),
                      (bounds[1][0], max(operator[1][1], bounds[1][1])))
        sites.append((part, bounds))
        x += GAP
    x0, x1 = min(xs), max(xs)
    ax.plot([0, x1], [0, 0], color="black", lw=1.0, zorder=0)
    # label height in glyph units
    label_height = (LABEL_PT / 72) * (x1 - x0) / width_in
    # drawn labels only
    y0 = min(ys + [y - label_height / 2 for y in label_ys])
    return sites, (x0, x1), (y0, max(ys))


def draw_arcs(ax, sites, glyph_top):
    """Regulatory arcs over the construct, stacked so no two overlap; returns the highest level."""
    # the lowest level that still leaves a descender below it on the tallest head
    base = glyph_top + CLEARANCE + HEAD_HEIGHT + DESCENT
    taken, top = {}, 0.0
    for sender, receiver, kind in RULES:
        senders = [index for index, (part, _) in enumerate(sites) if part == sender]
        if not senders:
            continue
        send_bounds = sites[senders[0]][1]
        send_x = (send_bounds[0][0] + send_bounds[1][0]) / 2
        # every copy of the site is regulated, not just the first
        for part, bounds in sites:
            if part != receiver:
                continue
            recv_x = (bounds[0][0] + bounds[1][0]) / 2
            span = (min(send_x, recv_x), max(send_x, recv_x))
            step = 0
            while not all(span[1] < start or span[0] > end
                          for start, end in taken.get(base + step * LEVEL_STEP, [])):
                step += 1
            level = base + step * LEVEL_STEP
            taken.setdefault(level, []).append(span)
            top = max(top, level)
            colour = COLOUR[sender]
            tip = bounds[1][1] + CLEARANCE
            foot = tip + HEAD_HEIGHT if kind == "stimulation" else tip
            ax.plot([send_x, send_x, recv_x, recv_x],
                    [send_bounds[1][1] + CLEARANCE, level, level, foot],
                    color=colour, lw=1.2, solid_joinstyle="miter", zorder=4)
            if kind == "stimulation":
                ax.add_patch(Polygon(
                    [(recv_x - HEAD_WIDTH / 2, foot), (recv_x + HEAD_WIDTH / 2, foot),
                     (recv_x, tip)],
                    closed=True, facecolor="white", edgecolor=colour, lw=1.2, zorder=4))
            else:
                ax.plot([recv_x - HEAD_WIDTH / 2, recv_x + HEAD_WIDTH / 2], [tip, tip],
                        color=colour, lw=1.8, zorder=4, solid_capstyle="butt")
    return top


def figure(parts, width_in):
    """One construct fitted to width_in, its height following so the glyphs stay undistorted."""
    fig = plt.figure(figsize=(width_in, 1.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    sites, (x0, x1), (y0, glyph_top) = draw(ax, parts, width_in)
    y1 = max(glyph_top, draw_arcs(ax, sites, glyph_top))
    x0, x1 = x0 - MARGIN, x1 + MARGIN
    y0, y1 = y0 - MARGIN, y1 + MARGIN_ABOVE
    fig.set_size_inches(width_in, width_in * (y1 - y0) / (x1 - x0))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    return fig


def composite(parts, chart_path, width_in, middle_path=None):
    """The construct drawn at width_in, stacked below the chart at chart_path.

    middle_path, when given, is laid between the two.
    """
    buffer = io.BytesIO()
    fig = figure(parts, width_in)
    fig.savefig(buffer, format="png", dpi=DPI)
    plt.close(fig)
    strip = Image.open(buffer)
    middle = Image.open(middle_path) if middle_path is not None else None
    with Image.open(chart_path) as chart:
        # widths must match; nothing is resampled
        assert strip.width == chart.width, (strip.width, chart.width)
        assert middle is None or middle.width == chart.width, (middle.width, chart.width)
        height = chart.height + strip.height + (middle.height if middle else 0)
        stacked = Image.new("RGB", (chart.width, height), "white")
        stacked.paste(chart, (0, 0))
        if middle is not None:
            stacked.paste(middle, (0, chart.height))
        stacked.paste(strip, (0, chart.height + (middle.height if middle else 0)))
    if middle is not None:
        middle.close()
    return stacked
