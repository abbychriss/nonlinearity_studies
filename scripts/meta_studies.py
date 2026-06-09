#!/usr/bin/env python3
"""
Config-driven meta-analysis of per-extension summary tables.

A "study" varies one independent variable (e.g. VR) across several
``extension_summary.csv`` tables produced by nonlinearity_studies and plots
chosen quantities against that variable, with one line per extension.

Datasets may optionally carry a ``series`` tag (a sub-category of the
independent variable, e.g. standard vs. increased deltaV at the same VR). Series
share the extension colour and are distinguished by linestyle; when there is
only one x value (each line is a lone point, so the linestyle is invisible) they
are distinguished by marker instead. They can be overlaid on one figure or split
onto separate figures.

Usage:
    python meta_studies.py -j config/VR_study.json

The X value of each table is supplied explicitly in the config (it is data, not
something parsed from the filepath), so the tool does not care how the source
files are named.
"""

import argparse
import json
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.size"] = 12

# Default colour per extension: the Tableau 10 palette.
DEFAULT_PALETTE = [
    "#4e79a7",  # blue
    "#f28e2b",  # orange
    "#e15759",  # red
    "#76b7b2",  # teal
    "#59a14f",  # green
    "#edc948",  # yellow
    "#b07aa1",  # purple
    "#ff9da7",  # pink
    "#9c755f",  # brown
    "#bab0ac",  # gray
]
DEFAULT_LINESTYLES = ["-", "--", "-.", ":"]
DEFAULT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

CONFIG_KEYS = {
    "study_name",
    "x_label",
    "invert_x",
    "output_dir",
    "quantities",
    "datasets",
    "series_label",
    "series_same_plot",
    "series_styles",
    "colors",
    "dpi",
    "show_plots",
    "fit_line",
    "save_report",
    "figsize",
    "plot_individual",
    "plot_together",
    "subplots_figsize",
    "subplots_ncols",
}

REQUIRED_KEYS = {"study_name", "x_label", "output_dir", "quantities", "datasets"}


def load_config(config_path):
    """Read and validate a study config, normalising dashes to underscores."""
    config_path = Path(config_path)
    with config_path.open("r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    if not isinstance(config, dict):
        raise ValueError(f"JSON config must contain an object at the top level: {config_path}")

    config = {key.replace("-", "_"): value for key, value in config.items()}

    unknown = sorted(set(config) - CONFIG_KEYS)
    if unknown:
        raise ValueError(
            f"Unknown config option(s): {', '.join(unknown)}. "
            f"Valid options: {', '.join(sorted(CONFIG_KEYS))}."
        )

    missing = sorted(REQUIRED_KEYS - set(config))
    if missing:
        raise ValueError(f"Missing required config option(s): {', '.join(missing)}.")

    # Defaults.
    config.setdefault("invert_x", False)
    config.setdefault("series_label", None)
    config.setdefault("series_same_plot", True)
    config.setdefault("series_styles", None)
    config.setdefault("colors", None)
    config.setdefault("dpi", 350)
    config.setdefault("show_plots", True)
    config.setdefault("fit_line", False)
    config.setdefault("save_report", True)
    config.setdefault("figsize", [9, 6.5])
    config.setdefault("plot_individual", True)
    config.setdefault("plot_together", False)
    config.setdefault("subplots_figsize", None)
    config.setdefault("subplots_ncols", None)

    for dataset in config["datasets"]:
        if "x" not in dataset or "table" not in dataset:
            raise ValueError(f"Each dataset needs an 'x' and a 'table': {dataset!r}")

    return config


def load_datasets(config):
    """Read every table, tag rows with x (and series), return one tidy frame."""
    records = []
    for dataset in config["datasets"]:
        df = pd.read_csv(dataset["table"]).sort_values("ext")
        df["x"] = dataset["x"]
        if "series" in dataset:
            df["series"] = dataset["series"]
        records.append(df)

    data = pd.concat(records, ignore_index=True)
    if "series" not in data.columns:
        data["series"] = None
    return data.sort_values(["ext", "series", "x"])


def _color_for(ext, extensions, colors):
    """Resolve an extension's colour from a list (cycled) or dict mapping."""
    if isinstance(colors, dict):
        # JSON keys are strings; accept both "1" and 1.
        return colors.get(str(ext), colors.get(ext, DEFAULT_PALETTE[0]))
    palette = colors if colors else DEFAULT_PALETTE
    return palette[extensions.index(ext) % len(palette)]


def _style_for(series, series_list, series_styles):
    """Resolve a series' linestyle from a dict mapping or a cycled default."""
    if series is None:
        return "-"
    if isinstance(series_styles, dict) and series in series_styles:
        return series_styles[series]
    cycle_styles = dict(zip(series_list, cycle(DEFAULT_LINESTYLES)))
    return cycle_styles[series]


def _marker_for(series, series_list, single_x):
    """Resolve a series' marker.

    Markers only differ between series when there is a single x value: each line
    is then a lone point and the linestyle is invisible, so the marker is the
    only thing distinguishing one series from another. Otherwise every line uses
    the default round marker and series are told apart by linestyle.
    """
    if not single_x or series is None:
        return "o"
    cycle_markers = dict(zip(series_list, cycle(DEFAULT_MARKERS)))
    return cycle_markers[series]


def _save_figure(config, column, suffix=""):
    name = f"{column}_vs_{config['study_name']}{suffix}.jpeg"
    out_path = Path(config["output_dir"]) / name
    plt.savefig(out_path, dpi=config["dpi"])
    return out_path


def _fit_enabled(spec, config):
    """Whether to fit a line for this quantity (per-quantity override of global)."""
    return spec.get("fit_line", config["fit_line"])


def compute_fits(data, config):
    """Fit y = slope*x + intercept for each (quantity, ext, series) line.

    Returns one record per line with the fit parameters, goodness of fit
    (Pearson R and R^2), and the percentage change in y from the first to the
    last x value (in the same reading order as the plot / printed table).

    Fitting is per quantity: a line is only fit when that quantity opts in via
    its spec's "fit_line" (falling back to the top-level "fit_line"). The
    percentage change is always recorded, since it doesn't depend on the fit.
    """
    extensions = sorted(data["ext"].unique())
    series_list = [s for s in data["series"].dropna().unique()] or [None]
    x_ascending = not config["invert_x"]
    records = []

    for column, spec in config["quantities"].items():
        fit_enabled = _fit_enabled(spec, config)
        for ext in extensions:
            for series in series_list:
                sub = data[data["ext"] == ext]
                if series is not None:
                    sub = sub[sub["series"] == series]
                sub = sub.dropna(subset=[column]).sort_values("x", ascending=x_ascending)
                if sub.empty:
                    continue

                xs = sub["x"].to_numpy(dtype=float)
                ys = sub[column].to_numpy(dtype=float)
                rec = {
                    "quantity": column, "ext": int(ext), "series": series,
                    "fit_enabled": fit_enabled,
                    "slope": None, "intercept": None, "r_value": None,
                    "r_squared": None, "pct_change": None,
                    "x_first": float(xs[0]), "x_last": float(xs[-1]),
                }
                if ys[0] != 0:
                    rec["pct_change"] = (ys[-1] - ys[0]) / abs(ys[0]) * 100.0
                if fit_enabled and len(np.unique(xs)) >= 2:
                    slope, intercept = np.polyfit(xs, ys, 1)
                    r = float(np.corrcoef(xs, ys)[0, 1])
                    rec.update(slope=float(slope), intercept=float(intercept),
                               r_value=r, r_squared=r * r)
                records.append(rec)

    return records


def _draw_line(ax, sub, column, color, style, label, fit_line, fit, marker="o"):
    """Draw one (ext, series) line onto ax: connected markers, or markers + fit."""
    xs = sub["x"].to_numpy(dtype=float)
    ys = sub[column].to_numpy(dtype=float)
    if fit_line:
        ax.plot(xs, ys, marker=marker, linestyle="none", color=color, label=label)
        if fit and fit["slope"] is not None:
            x_ends = np.array([xs.min(), xs.max()])
            fit_label = f"{label} fit (R² = {fit['r_squared']:.3f})"
            ax.plot(x_ends, fit["slope"] * x_ends + fit["intercept"],
                    color=color, linestyle=style, alpha=0.5, label=fit_label)
    else:
        ax.plot(xs, ys, marker=marker, color=color, linestyle=style, label=label)


def _render_quantity(ax, subset, column, spec, config, fits, series_list,
                     group_series, legend_fontsize=None):
    """Draw one quantity onto an axes.

    ``subset`` is the data to plot (the full set, or one series in split mode);
    ``group_series`` names that single series when the axes holds just one, else
    None (which, with multiple series, means overlay them all on this axes).
    """
    extensions = sorted(subset["ext"].unique())
    has_series = len(series_list) > 0
    multi_series = len(series_list) > 1
    overlay = multi_series and group_series is None
    fit_line = _fit_enabled(spec, config)
    fit_lookup = {(f["ext"], f["series"]): f for f in fits if f["quantity"] == column}
    # A single x value makes the linestyle invisible (each line is a lone point),
    # so series are distinguished by marker instead.
    single_x = subset["x"].nunique() <= 1

    for ext in extensions:
        ext_rows = subset[subset["ext"] == ext]
        color = _color_for(ext, extensions, config["colors"])
        if overlay:
            for series in series_list:
                sub = ext_rows[ext_rows["series"] == series].sort_values("x")
                if sub.empty:
                    continue
                style = _style_for(series, series_list, config["series_styles"])
                marker = _marker_for(series, series_list, single_x)
                _draw_line(ax, sub, column, color, style, f"EXT{ext} ({series})",
                           fit_line, fit_lookup.get((ext, series)), marker=marker)
        else:
            sub = ext_rows.sort_values("x")
            if sub.empty:
                continue
            series = group_series if group_series is not None else (
                series_list[0] if has_series else None)
            style = _style_for(series, series_list, config["series_styles"])
            marker = _marker_for(series, series_list, single_x)
            _draw_line(ax, sub, column, color, style, f"EXT{ext}",
                       fit_line, fit_lookup.get((ext, series)), marker=marker)

    # With a single series on the axes, name it in the title rather than
    # repeating it on every legend entry (e.g. "Gain vs VR (deltaV = Standard)").
    title = spec.get("title", f"{column} vs {config['study_name']}")
    if not overlay:
        fig_series = group_series if group_series is not None else (
            series_list[0] if has_series else None)
        if fig_series is not None:
            title = f"{title} ({config['series_label'] or 'series'} = {fig_series})"

    if config["invert_x"]:
        ax.invert_xaxis()
    ax.set_xlabel(config["x_label"])
    ax.set_ylabel(spec.get("ylabel", column))
    ax.set_title(title)
    ax.legend(fontsize=legend_fontsize)
    ax.grid(True, alpha=0.3)


def plot_quantity(data, column, spec, config, fits):
    """Plot one quantity as its own figure(s). Returns saved figure paths."""
    series_list = [s for s in data["series"].dropna().unique()]
    multi_series = len(series_list) > 1
    saved = []

    # Split onto separate figures only when several series can't share one plot.
    if multi_series and not config["series_same_plot"]:
        groups = [(s, data[data["series"] == s], f"_{s}") for s in series_list]
    else:
        groups = [(None, data, "")]

    for group_series, subset, file_suffix in groups:
        fig, ax = plt.subplots(figsize=config["figsize"])
        _render_quantity(ax, subset, column, spec, config, fits, series_list, group_series)
        fig.tight_layout()
        saved.append(_save_figure(config, column, file_suffix))

    return saved


def plot_all_together(data, config, fits):
    """Plot every quantity as subplots in one figure. Returns the saved path."""
    series_list = [s for s in data["series"].dropna().unique()]
    quantities = list(config["quantities"].items())
    n = len(quantities)
    # Default to a roughly square grid so cells stay large (e.g. 3 -> 2x2).
    ncols = config["subplots_ncols"] or int(np.ceil(np.sqrt(n)))
    nrows = -(-n // ncols)  # ceil
    figsize = config["subplots_figsize"] or [6.5 * ncols, 5.5 * nrows]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    flat_axes = axes.flatten()
    for ax, (column, spec) in zip(flat_axes, quantities):
        # Subplots overlay series in each cell (one cell per quantity); shrink the
        # legend so it doesn't crowd out the data at subplot scale.
        _render_quantity(ax, data, column, spec, config, fits, series_list, None,
                         legend_fontsize="small")
    for ax in flat_axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    out_path = Path(config["output_dir"]) / f"{config['study_name']}_subplots.jpeg"
    fig.savefig(out_path, dpi=config["dpi"])
    return out_path


def write_report(config, table_str, fits):
    """Write the summary table plus per-line fit stats to a text report."""
    lines = [
        f"Study: {config['study_name']}",
        f"x-axis: {config['x_label']}",
        "",
        "=== Summary table ===",
        table_str,
        "",
        "=== Line fits (y = slope * x + intercept) ===",
    ]
    for column, spec in config["quantities"].items():
        lines.append("")
        lines.append(f"{spec.get('ylabel', column)}  [{column}]")
        for fit in (r for r in fits if r["quantity"] == column):
            tag = f"EXT{fit['ext']}"
            if fit["series"] is not None:
                tag += f" ({fit['series']})"
            if fit["pct_change"] is not None:
                delta = (f"%change(x={fit['x_first']:g} -> {fit['x_last']:g}) = "
                         f"{fit['pct_change']:+.2f}%")
            else:
                delta = "%change = n/a"
            if not fit["fit_enabled"]:
                parts = ["fit not requested", delta]
            elif fit["slope"] is not None:
                parts = [f"y = {fit['slope']:.4g} x {fit['intercept']:+.4g}",
                         f"R = {fit['r_value']:.4f}, R^2 = {fit['r_squared']:.4f}",
                         delta]
            else:
                parts = ["fit unavailable (need >=2 distinct x)", "R = n/a", delta]
            lines.append(f"  {tag}: {' | '.join(parts)}")

    report_path = Path(config["output_dir"]) / f"{config['study_name']}_report.txt"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def write_fit_stats_csv(config, fits):
    """Write the per-line fit stats as a machine-readable CSV."""
    columns = ["quantity", "ext", "series", "fit_enabled", "slope", "intercept",
               "r_value", "r_squared", "pct_change", "x_first", "x_last"]
    csv_path = Path(config["output_dir"]) / f"{config['study_name']}_fit_stats.csv"
    pd.DataFrame(fits, columns=columns).to_csv(csv_path, index=False)
    return csv_path


def run_study(config):
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)
    data = load_datasets(config)

    # Order the printed table the same way the plots read: x descending when
    # the x-axis is inverted (e.g. high VR to low VR), ascending otherwise.
    x_ascending = not config["invert_x"]
    sort_cols = ["series", "x", "ext"] if config["series_label"] else ["x", "ext"]
    ascending = [True if col != "x" else x_ascending for col in sort_cols]
    # Display the x and series columns under their human-readable labels.
    renames = {"x": config["x_label"]}
    if config["series_label"]:
        renames["series"] = config["series_label"]
    table_str = (data.sort_values(sort_cols, ascending=ascending)
                 .rename(columns=renames).to_string(index=False))
    print(table_str)

    fits = compute_fits(data, config)

    if config["plot_individual"]:
        for column, spec in config["quantities"].items():
            for path in plot_quantity(data, column, spec, config, fits):
                print(f"Saved {path}")

    if config["plot_together"]:
        print(f"Saved {plot_all_together(data, config, fits)}")

    if config["save_report"]:
        print(f"Saved {write_report(config, table_str, fits)}")
        print(f"Saved {write_fit_stats_csv(config, fits)}")

    if config["show_plots"]:
        plt.show()


def init_argparse(args=None):
    parser = argparse.ArgumentParser(
        description="Plot per-extension summary quantities across a varied parameter."
    )
    parser.add_argument(
        "-j", "--json", dest="config", required=True,
        help="Path to the study JSON config.",
    )
    return parser.parse_args(args)


def main(args=None):
    parsed = init_argparse(args)
    config = load_config(parsed.config)
    run_study(config)


if __name__ == "__main__":
    main()
