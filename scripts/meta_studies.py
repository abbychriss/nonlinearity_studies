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
onto separate figures (``series_same_plot``). When split, set
``"uniform_series_style": true`` so every per-series figure uses the same solid
line and circle marker instead of cycling styles that would make the separate
figures look needlessly different.

Set ``"x_axis": "extension"`` to flip the orientation: extension goes on the
x-axis, each series becomes its own coloured line, and the independent-variable
value (e.g. VR = -6) is pinned in the title. One figure is produced per distinct
x value. This is the natural view when the independent variable is held fixed
and you want to compare a quantity across extensions for a few series.

Usage:
    python meta_studies.py -j config/VR_study.json

The X value of each table is supplied explicitly in the config (it is data, not
something parsed from the filepath), so the tool does not care how the source
files are named.
"""

import argparse
import json
import re
from glob import glob
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.size"] = 12

# Default colour per extension: the Tableau 10 palette.
DEFAULT_PALETTE = [
    "#4e79a7",  # blue
    "#e15759",  # red
    "#76b7b2",  # teal
    "#f28e2b",  # orange
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
    "x_axis",
    "invert_x",
    "output_dir",
    "quantities",
    "datasets",
    "series_label",
    "series_column",
    "series_same_plot",
    "series_line_styles",
    "series_markers",
    "uniform_series_style",
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
    config.setdefault("x_axis", "value")
    config.setdefault("invert_x", False)
    config.setdefault("series_label", None)
    config.setdefault("series_column", None)
    config.setdefault("series_same_plot", True)
    config.setdefault("series_line_styles", None)
    config.setdefault("series_markers", None)
    config.setdefault("uniform_series_style", False)
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

    if config["x_axis"] not in ("value", "extension"):
        raise ValueError(
            f"'x_axis' must be 'value' or 'extension', got {config['x_axis']!r}."
        )

    for dataset in config["datasets"]:
        if "x" not in dataset or "table" not in dataset:
            raise ValueError(f"Each dataset needs an 'x' and a 'table': {dataset!r}")
        if "exclude_ext" in dataset and not isinstance(dataset["exclude_ext"], list):
            raise ValueError(
                f"'exclude_ext' must be a list of extension numbers: {dataset!r}"
            )

    return config


def _resolve_quantity_columns(df, quantities, table):
    """Ensure each configured quantity is present as a column named by its key.

    By default a quantity key must match a CSV header. A quantity may instead
    name its source column explicitly via an absolute 0-based ``"column"`` index
    into the CSV; the addressed column's values are then exposed under the
    quantity key so the rest of the pipeline can refer to it by key regardless of
    the header's actual name. This makes it easy to pick out just the columns you
    want to plot without retyping their headers.
    """
    n_cols = len(df.columns)
    for key, spec in quantities.items():
        index = spec.get("column") if isinstance(spec, dict) else None
        if index is None:
            if key not in df.columns:
                raise ValueError(
                    f"Quantity {key!r} is not a column in {table} and no 'column' "
                    f"index was given. Available columns: {', '.join(map(str, df.columns))}."
                )
            continue
        if not isinstance(index, int) or isinstance(index, bool):
            raise ValueError(
                f"Quantity {key!r} 'column' must be an integer index, got {index!r}."
            )
        if not -n_cols <= index < n_cols:
            raise ValueError(
                f"Quantity {key!r} 'column' index {index} is out of range for {table} "
                f"({n_cols} columns: {', '.join(map(str, df.columns))})."
            )
        df[key] = df.iloc[:, index]
    return df


def _resolve_table_path(table):
    """Expand a dataset table path, which may contain glob wildcards (``*``, ``?``).

    The path is supplied by the config and often points at an auto-named plots
    subdirectory (e.g. ``plots/avg*/extension_summary.csv``), so a wildcard saves
    the user from pasting the full hashed directory name. Exactly one file must
    match: zero matches or an ambiguous several are raised as errors rather than
    silently picking one. A path with no wildcard is returned as-is when it exists.
    """
    matches = sorted(glob(table, recursive=True))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(f"No file matches table path: {table!r}")
    raise ValueError(
        f"Table path {table!r} is ambiguous; {len(matches)} files match: "
        f"{', '.join(matches)}. Make the pattern more specific."
    )


def load_datasets(config):
    """Read every table, tag rows with x (and series), return one tidy frame.

    A dataset may carry ``"exclude_ext": [...]`` to drop specific extensions from
    that dataset only (e.g. an extension that misbehaves at one operating point),
    leaving those extensions intact in the other datasets.

    The series dimension comes from one of two places. Normally it is a per-dataset
    tag (``"series": "-8"``), constant across the whole table. Alternatively, when
    ``"series_column"`` is set in the config, the series is read from that CSV
    column instead, so a table that varies a second parameter within itself (e.g.
    ``resolution_summary.csv`` with one row per charge ``q``) is split into one
    series per distinct value of that column. The two are mutually exclusive;
    ``series_column`` takes precedence if both are present.
    """
    series_column = config.get("series_column")
    records = []
    for dataset in config["datasets"]:
        table = _resolve_table_path(dataset["table"])
        df = pd.read_csv(table).sort_values("ext")
        df = _resolve_quantity_columns(df, config["quantities"], table)
        exclude_ext = dataset.get("exclude_ext")
        if exclude_ext:
            df = df[~df["ext"].isin(exclude_ext)]
        df["x"] = dataset["x"]
        if series_column:
            if series_column not in df.columns:
                raise ValueError(
                    f"series_column {series_column!r} not found in {dataset['table']}. "
                    f"Available columns: {', '.join(map(str, df.columns))}."
                )
            df["series"] = df[series_column]
        elif "series" in dataset:
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


def _style_for(series, series_list, series_line_styles, uniform=False):
    """Resolve a series' linestyle from a dict mapping or a cycled default.

    With ``uniform`` set, every series gets the default solid line; this is for
    split plots (one series per figure), where cycling styles only makes the
    separate figures look gratuitously different.
    """
    if series is None or uniform:
        return "-"
    if isinstance(series_line_styles, dict) and series in series_line_styles:
        return series_line_styles[series]
    cycle_styles = dict(zip(series_list, cycle(DEFAULT_LINESTYLES)))
    return cycle_styles[series]


def _marker_for(series, series_list, single_x, series_markers=None, uniform=False):
    """Resolve a series' marker from a dict mapping or a default.

    An explicit ``series_markers`` mapping always wins: a series named there gets
    its configured marker on every figure. Otherwise markers only differ between
    series when there is a single x value (each line is then a lone point and the
    linestyle is invisible, so the marker is the only distinguisher); with several
    x values every line uses the default round marker and series are told apart by
    linestyle. With ``uniform`` set, every series gets the default marker.
    """
    if uniform or series is None:
        return "o"
    if isinstance(series_markers, dict) and series in series_markers:
        return series_markers[series]
    if not single_x:
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
                style = _style_for(series, series_list, config["series_line_styles"])
                marker = _marker_for(series, series_list, single_x,
                                     config["series_markers"])
                tag = _series_tag(config["series_label"] or "series", series)
                _draw_line(ax, sub, column, color, style, f"EXT{ext} ({tag})",
                           fit_line, fit_lookup.get((ext, series)), marker=marker)
        else:
            sub = ext_rows.sort_values("x")
            if sub.empty:
                continue
            series = group_series if group_series is not None else (
                series_list[0] if has_series else None)
            uniform = config["uniform_series_style"]
            style = _style_for(series, series_list, config["series_line_styles"],
                               uniform=uniform)
            marker = _marker_for(series, series_list, single_x,
                                 config["series_markers"], uniform=uniform)
            _draw_line(ax, sub, column, color, style, f"EXT{ext}",
                       fit_line, fit_lookup.get((ext, series)), marker=marker)

    # With a single series on the axes, name it in the title rather than
    # repeating it on every legend entry (e.g. "Gain vs VR (deltaV = Standard)").
    title = spec.get("title", f"{column} vs {config['study_name']}")
    if not overlay:
        fig_series = group_series if group_series is not None else (
            series_list[0] if has_series else None)
        if fig_series is not None:
            title = f"{title} ({_labeled_value(config['series_label'] or 'series', fig_series)})"

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


def _series_tag(label, value):
    """Compact series value for a legend entry.

    Numeric series are formatted without trailing zeros and, when the label
    carries a unit (``"Charge (e-)"``), gain that unit directly (200 -> "200e-");
    the label name itself is omitted since the legend already groups by series.
    Non-numeric series (e.g. ``"Standard"``) are used verbatim.
    """
    try:
        text = f"{float(value):g}"
    except (TypeError, ValueError):
        return str(value)
    if label:
        match = re.fullmatch(r"\s*(.*?)\s*\(([^)]*)\)\s*", label)
        if match:
            text += match.group(2)
    return text


def _labeled_value(label, value):
    """Render ``"<name> = <value><unit>"`` from a label like ``"VR (V)"``.

    A label of the form ``"NAME (UNIT)"`` attaches the unit directly to the value
    (e.g. ``"VR (V)"`` with -8 -> ``"VR = -8V"``); otherwise the label is used
    verbatim (``"<label> = <value>"``). Numeric values (including numeric strings
    like ``"-8"`` from a JSON series tag) are formatted compactly; non-numeric
    values (e.g. ``"Standard"``) are used as-is.
    """
    try:
        value = f"{float(value):g}"
    except (TypeError, ValueError):
        value = str(value)
    match = re.fullmatch(r"\s*(.*?)\s*\(([^)]*)\)\s*", label)
    if match:
        name, unit = match.group(1), match.group(2)
        return f"{name} = {value}{unit}"
    return f"{label} = {value}"


def _render_quantity_vs_ext(ax, subset, column, spec, config, x_value,
                            series_list, legend_fontsize=None):
    """Draw one quantity with extension on the x-axis at a fixed x value.

    Here the roles are flipped relative to ``_render_quantity``: extension is the
    independent axis and each series is its own colour and marker, while the
    held-fixed independent-variable value (e.g. VR = -6) is named in the title.
    Extensions are independent channels, not a continuum, so points are drawn
    unconnected (no line between extensions).
    """
    series_to_plot = series_list if series_list else [None]
    palette = (config["colors"] if isinstance(config["colors"], list)
               and config["colors"] else DEFAULT_PALETTE)

    for idx, series in enumerate(series_to_plot):
        sub = subset if series is None else subset[subset["series"] == series]
        sub = sub.dropna(subset=[column]).sort_values("ext")
        if sub.empty:
            continue
        color = palette[idx % len(palette)]
        marker = DEFAULT_MARKERS[idx % len(DEFAULT_MARKERS)] if series is not None else "o"
        label = str(series) if series is not None else spec.get("ylabel", column)
        ax.plot(sub["ext"].to_numpy(dtype=float), sub[column].to_numpy(dtype=float),
                marker=marker, linestyle="none", color=color, label=label)

    base_title = spec.get("title_vs_ext", f"{spec.get('ylabel', column)} vs extension")
    ax.set_title(f"{base_title} ({_labeled_value(config['x_label'], x_value)})")
    ax.set_xlabel("Extension")
    ax.set_ylabel(spec.get("ylabel", column))
    ax.set_xticks(sorted(subset["ext"].unique()))
    ax.legend(fontsize=legend_fontsize)
    ax.grid(True, alpha=0.3)


def plot_quantity_vs_ext(data, column, spec, config):
    """Plot one quantity vs extension, one figure per distinct x value.

    Returns the saved figure paths.
    """
    series_list = [s for s in data["series"].dropna().unique()]
    multi_x = data["x"].nunique() > 1
    saved = []
    for x_value in sorted(data["x"].unique()):
        subset = data[data["x"] == x_value]
        fig, ax = plt.subplots(figsize=config["figsize"])
        _render_quantity_vs_ext(ax, subset, column, spec, config, x_value, series_list)
        fig.tight_layout()
        suffix = f"_x{x_value:g}" if multi_x else ""
        out_path = (Path(config["output_dir"])
                    / f"{column}_vs_ext_{config['study_name']}{suffix}.jpeg")
        fig.savefig(out_path, dpi=config["dpi"])
        saved.append(out_path)
    return saved


def plot_all_together_vs_ext(data, config):
    """Plot every quantity vs extension as subplots, one figure per x value.

    Returns the saved figure paths.
    """
    series_list = [s for s in data["series"].dropna().unique()]
    quantities = list(config["quantities"].items())
    n = len(quantities)
    ncols = config["subplots_ncols"] or int(np.ceil(np.sqrt(n)))
    nrows = -(-n // ncols)  # ceil
    figsize = config["subplots_figsize"] or [6.5 * ncols, 5.5 * nrows]
    multi_x = data["x"].nunique() > 1
    saved = []

    for x_value in sorted(data["x"].unique()):
        subset = data[data["x"] == x_value]
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        flat_axes = axes.flatten()
        for ax, (column, spec) in zip(flat_axes, quantities):
            _render_quantity_vs_ext(ax, subset, column, spec, config, x_value,
                                    series_list, legend_fontsize="small")
        for ax in flat_axes[n:]:
            ax.set_visible(False)
        fig.tight_layout()
        suffix = f"_x{x_value:g}" if multi_x else ""
        out_path = (Path(config["output_dir"])
                    / f"{config['study_name']}_ext_subplots{suffix}.jpeg")
        fig.savefig(out_path, dpi=config["dpi"])
        saved.append(out_path)
    return saved


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
    vs_extension = config["x_axis"] == "extension"

    if config["plot_individual"]:
        for column, spec in config["quantities"].items():
            paths = (plot_quantity_vs_ext(data, column, spec, config) if vs_extension
                     else plot_quantity(data, column, spec, config, fits))
            for path in paths:
                print(f"Saved {path}")

    if config["plot_together"]:
        paths = (plot_all_together_vs_ext(data, config) if vs_extension
                 else [plot_all_together(data, config, fits)])
        for path in paths:
            print(f"Saved {path}")

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
