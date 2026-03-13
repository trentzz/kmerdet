#!/usr/bin/env python3
"""plot_results.py — Generate figures from kmerdet benchmark metrics.

Produces the following figures:
  1. sensitivity_by_vaf.png      — Sensitivity vs VAF curve (ROC-like)
  2. runtime_scaling.png         — Runtime vs number of targets (linear + log)
  3. memory_by_dataset.png       — Peak RSS memory by dataset size
  4. variant_type_breakdown.png  — Stacked bar chart by variant type
  5. multi_k_comparison.png      — Multi-k vs single-k sensitivity comparison
  6. filtering_venn.png          — Before/after filtering Venn diagrams
  7. sensitivity_heatmap.png     — Heatmap: variant type × VAF range
  8. threshold_sweep.png         — Sensitivity/precision vs rVAF threshold

Usage:
    python3 plot_results.py --metrics analysis_output/metrics.json
    python3 plot_results.py --results-dir /path/to/benchmark/results
    python3 plot_results.py --per-type per_type.tsv --per-vaf per_vaf_bin.tsv
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from typing import Optional

import matplotlib
matplotlib.use("Agg")   # non-interactive backend; safe in CI
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------

THESIS_COLOR     = "#e74c3c"    # red — km/thesis baseline
KMERDET_COLOR    = "#2980b9"    # blue — kmerdet
NEUTRAL_COLOR    = "#95a5a6"    # grey — supplementary info

VARIANT_COLORS = {
    "Substitution": "#3498db",
    "Insertion":    "#2ecc71",
    "Deletion":     "#e74c3c",
    "ITD":          "#f39c12",
    "Complex":      "#9b59b6",
}

VAF_BIN_LABELS = ["0–0.1%", "0.1–1%", "1–5%", "5–100%"]

# Thesis baseline sensitivity values for annotation
THESIS_BASELINES = {
    "Substitution": 0.77,
    "Insertion":    0.38,
    "Deletion":     0.38,
}

FIG_DPI = 150
FIG_SIZE_SINGLE = (7, 5)
FIG_SIZE_WIDE   = (10, 5)
FIG_SIZE_SQUARE = (6, 6)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(v) -> Optional[float]:
    if v is None:
        return None
    try:
        f = float(v)
        return None if math.isnan(f) else f
    except (TypeError, ValueError):
        return None


def savefig(fig: plt.Figure, path: pathlib.Path) -> None:
    fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 1. Sensitivity vs VAF curve
# ---------------------------------------------------------------------------

def plot_sensitivity_by_vaf(
    per_vaf_df: pd.DataFrame,
    out_path: pathlib.Path,
    label: str = "kmerdet",
) -> None:
    """Bar chart of sensitivity at each VAF bin, with thesis baseline overlay."""
    if per_vaf_df.empty:
        print("  Skipping: no VAF bin data.")
        return

    bins = per_vaf_df[per_vaf_df["present"] > 0].copy()
    if bins.empty:
        print("  Skipping: no non-empty VAF bins.")
        return

    n = len(bins)
    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=FIG_SIZE_WIDE)

    senss = [_safe_float(s) or 0.0 for s in bins["sensitivity"]]
    bars = ax.bar(x, senss, width, label=label, color=KMERDET_COLOR, alpha=0.85)

    # Annotate bars
    for bar, s in zip(bars, senss):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.01,
            f"{s:.2f}",
            ha="center", va="bottom", fontsize=9,
        )

    # Thesis baseline horizontal lines (if we have per-bin baselines or use SNV/INDEL overall)
    thesis_snv_line  = ax.axhline(0.77, ls="--", color=THESIS_COLOR, lw=1.2,
                                   label="km SNV baseline (77%)", alpha=0.7)
    thesis_indel_line = ax.axhline(0.38, ls=":", color=THESIS_COLOR, lw=1.2,
                                    label="km INDEL baseline (38%)", alpha=0.7)

    bin_labels = list(bins.get("bin_label", bins.index))
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, fontsize=10)
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("True VAF bin", fontsize=11)
    ax.set_ylabel("Sensitivity (TP / (TP + FN))", fontsize=11)
    ax.set_title("Detection Sensitivity by Variant Allele Frequency", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")

    # Annotate count above x-axis
    for i, (_, row) in enumerate(bins.iterrows()):
        ax.text(x[i], -0.07, f"n={int(row['present'])}",
                ha="center", va="top", fontsize=8, color="grey",
                transform=ax.get_xaxis_transform())

    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 2. Runtime scaling (targets vs time)
# ---------------------------------------------------------------------------

def plot_runtime_scaling(
    timing_records: list[dict],
    out_path: pathlib.Path,
) -> None:
    """Two-panel plot: runtime vs targets on linear and log scales."""
    if not timing_records:
        print("  Skipping: no timing data.")
        return

    df = pd.DataFrame(timing_records)
    df = df.dropna(subset=["wall_sec"])

    if df.empty:
        print("  Skipping: all timing records are null.")
        return

    fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=FIG_SIZE_WIDE, sharey=False)

    for ax in (ax_lin, ax_log):
        ax.scatter(df.index, df["wall_sec"], color=KMERDET_COLOR, zorder=3)
        if len(df) > 1:
            ax.plot(df.index, df["wall_sec"], color=KMERDET_COLOR, alpha=0.6)

        # Thesis reference line at 120 s (2 min detect step for 50 targets)
        ax.axhline(120, ls="--", color=THESIS_COLOR, lw=1.0,
                   label="km detect ~120 s (thesis)", alpha=0.7)

        ax.set_xlabel("Benchmark run", fontsize=10)
        ax.set_ylabel("Wall-clock time (s)", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    ax_log.set_yscale("log")
    ax_log.set_ylabel("Wall-clock time (s, log scale)", fontsize=10)
    ax_lin.set_title("Runtime (linear)", fontsize=11)
    ax_log.set_title("Runtime (log scale)", fontsize=11)

    fig.suptitle("kmerdet Detection Runtime", fontsize=12)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 3. Memory usage by dataset
# ---------------------------------------------------------------------------

def plot_memory_usage(
    timing_records: list[dict],
    out_path: pathlib.Path,
) -> None:
    """Bar chart of peak RSS memory by benchmark run."""
    df = pd.DataFrame(timing_records)
    df = df.dropna(subset=["peak_rss_mb"])

    if df.empty:
        print("  Skipping: no memory data.")
        return

    fig, ax = plt.subplots(figsize=FIG_SIZE_SINGLE)
    x = np.arange(len(df))
    bars = ax.bar(x, df["peak_rss_mb"], color=KMERDET_COLOR, alpha=0.85)

    for bar, v in zip(bars, df["peak_rss_mb"]):
        ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + 5,
                f"{v:.0f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(
        df.get("label", df.index).tolist(),
        rotation=30, ha="right", fontsize=9,
    )
    ax.set_ylabel("Peak RSS (MB)", fontsize=11)
    ax.set_title("Peak Memory Usage by Benchmark Run", fontsize=12)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 4. Variant type breakdown (stacked bar)
# ---------------------------------------------------------------------------

def plot_variant_type_breakdown(
    per_type_df: pd.DataFrame,
    out_path: pathlib.Path,
) -> None:
    """Stacked bar: TP (detected) vs FN (missed) per variant type."""
    if per_type_df.empty:
        print("  Skipping: no per-type data.")
        return

    df = per_type_df[per_type_df["present"] > 0].copy()
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=FIG_SIZE_WIDE)
    x = np.arange(len(df))
    width = 0.5

    tp_vals = df["tp"].astype(int)
    fn_vals = df["fn"].astype(int) if "fn" in df.columns else (
        df["present"].astype(int) - tp_vals
    )

    ax.bar(x, tp_vals, width, label="Detected (TP)", color=KMERDET_COLOR, alpha=0.85)
    ax.bar(x, fn_vals, width, bottom=tp_vals, label="Missed (FN)", color=THESIS_COLOR, alpha=0.75)

    # Overlay sensitivity text
    for i, (_, row) in enumerate(df.iterrows()):
        s = _safe_float(row["sensitivity"])
        if s is not None:
            ax.text(x[i], int(row["present"]) + 0.5, f"{s:.0%}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

        # Thesis baseline marker
        vtype = row["variant_type"]
        if vtype in THESIS_BASELINES:
            baseline_n = THESIS_BASELINES[vtype] * int(row["present"])
            ax.plot([x[i] - 0.3, x[i] + 0.3], [baseline_n, baseline_n],
                    ls="--", color=THESIS_COLOR, lw=2.0)

    ax.set_xticks(x)
    ax.set_xticklabels(df["variant_type"].tolist(), fontsize=10)
    ax.set_ylabel("Number of variants", fontsize=11)
    ax.set_title("Detection Performance by Variant Type", fontsize=12)
    ax.legend(fontsize=9)

    # Add thesis baseline legend entry
    thesis_patch = mpatches.Patch(color=THESIS_COLOR, alpha=0.0, label="")
    ax.add_patch(thesis_patch)
    ax.plot([], [], ls="--", color=THESIS_COLOR, lw=2.0, label="km baseline sensitivity")
    ax.legend(fontsize=9)

    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 5. Multi-k vs single-k comparison
# ---------------------------------------------------------------------------

def plot_multi_k_comparison(
    single_k_type: pd.DataFrame,
    multi_k_type: pd.DataFrame,
    out_path: pathlib.Path,
) -> None:
    """Grouped bar chart comparing single-k vs multi-k sensitivity by type."""
    if single_k_type.empty and multi_k_type.empty:
        print("  Skipping: no multi-k data.")
        return

    # Merge on variant_type
    merged = single_k_type[["variant_type", "sensitivity"]].rename(
        columns={"sensitivity": "single_k"}
    ).merge(
        multi_k_type[["variant_type", "sensitivity"]].rename(
            columns={"sensitivity": "multi_k"}
        ),
        on="variant_type",
        how="outer",
    )
    merged = merged.fillna(0)

    n = len(merged)
    x = np.arange(n)
    width = 0.3

    fig, ax = plt.subplots(figsize=FIG_SIZE_WIDE)

    bars1 = ax.bar(x - width/2, merged["single_k"], width,
                   label="Single-k (k=31)", color=NEUTRAL_COLOR, alpha=0.85)
    bars2 = ax.bar(x + width/2, merged["multi_k"], width,
                   label="Multi-k (k=21+31)", color=KMERDET_COLOR, alpha=0.85)

    def annotate(bars, vals):
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + 0.01,
                        f"{v:.2f}", ha="center", va="bottom", fontsize=8)

    annotate(bars1, merged["single_k"])
    annotate(bars2, merged["multi_k"])

    ax.set_xticks(x)
    ax.set_xticklabels(merged["variant_type"].tolist(), fontsize=10)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Sensitivity", fontsize=11)
    ax.set_title("Multi-k vs Single-k Detection Sensitivity", fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 6. Before/after filtering Venn diagram (approximate)
# ---------------------------------------------------------------------------

def plot_filtering_venn(
    n_before: int,
    n_after: int,
    n_tp: int,
    n_fp_removed: int,
    n_fn_added: int,
    out_path: pathlib.Path,
) -> None:
    """Simple area-proportional circles showing filter effect."""
    if n_before == 0:
        print("  Skipping: no filtering data.")
        return

    fig, ax = plt.subplots(figsize=FIG_SIZE_SQUARE)
    ax.set_aspect("equal")
    ax.set_xlim(-2, 4)
    ax.set_ylim(-2, 2)
    ax.axis("off")

    scale = math.sqrt(n_before / 200) if n_before > 0 else 1.0
    r_before = 1.2 * scale
    r_after  = max(0.3, r_before * math.sqrt(n_after / max(n_before, 1)))

    # Before circle
    circle_before = plt.Circle((0, 0), r_before, color=NEUTRAL_COLOR, alpha=0.35, lw=2)
    ax.add_patch(circle_before)
    ax.text(0, r_before + 0.15, f"Before filtering\n(n={n_before})",
            ha="center", va="bottom", fontsize=10)

    # After circle (overlapping)
    overlap_offset = r_before * 0.4
    circle_after = plt.Circle((overlap_offset, 0), r_after, color=KMERDET_COLOR, alpha=0.45, lw=2)
    ax.add_patch(circle_after)
    ax.text(overlap_offset + r_after + 0.1, 0,
            f"After filtering\n(n={n_after})", ha="left", va="center", fontsize=10)

    # Annotations
    ax.text(overlap_offset, -r_after - 0.2,
            f"TP kept: {n_tp}\nFP removed: {n_fp_removed}\nFN added: {n_fn_added}",
            ha="center", va="top", fontsize=9, color="grey")

    ax.set_title("Effect of Filtering on Variant Calls", fontsize=12, pad=20)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 7. Sensitivity heatmap: variant type × VAF bin
# ---------------------------------------------------------------------------

def plot_sensitivity_heatmap(
    per_variant_df: pd.DataFrame,
    out_path: pathlib.Path,
    vaf_bins: Optional[list[tuple[float, float, str]]] = None,
) -> None:
    """Heatmap of sensitivity: rows = variant type, columns = VAF bin."""
    if per_variant_df.empty:
        print("  Skipping: no per-variant data.")
        return

    if vaf_bins is None:
        vaf_bins = [
            (0.0,    0.001,  "0–0.1%"),
            (0.001,  0.01,   "0.1–1%"),
            (0.01,   0.05,   "1–5%"),
            (0.05,   1.0,    "5–100%"),
        ]

    present = per_variant_df[per_variant_df.get("true_vaf", pd.Series()).gt(0)].copy()
    if "true_vaf" not in present.columns or present.empty:
        print("  Skipping: no present variants in per-variant data.")
        return

    present["true_vaf"] = pd.to_numeric(present["true_vaf"], errors="coerce")
    present = present.dropna(subset=["true_vaf", "variant_type"])
    present = present[present["true_vaf"] > 0]

    all_types = [t for t in VARIANT_TYPE_ORDER if t in present["variant_type"].unique()]
    all_types += [t for t in present["variant_type"].unique() if t not in VARIANT_TYPE_ORDER]

    # Build sensitivity matrix
    matrix = np.full((len(all_types), len(vaf_bins)), float("nan"))

    for j, (low, high, _) in enumerate(vaf_bins):
        in_bin = present[
            (present["true_vaf"] >= low) &
            (present["true_vaf"] < high if high < 1.0 else present["true_vaf"] <= high)
        ]
        for i, vtype in enumerate(all_types):
            sub = in_bin[in_bin["variant_type"] == vtype]
            if len(sub) == 0:
                continue
            n_tp = int((sub["status"] == "TP").sum())
            matrix[i, j] = n_tp / len(sub)

    bin_labels = [b[2] for b in vaf_bins]
    fig, ax = plt.subplots(figsize=(max(6, len(vaf_bins) * 1.5), max(4, len(all_types) * 0.8)))

    # Replace NaN with -1 for plotting (shown as grey)
    plot_matrix = np.where(np.isnan(matrix), -1, matrix)
    cmap = matplotlib.cm.get_cmap("RdYlGn").copy()
    cmap.set_under("lightgrey")

    im = ax.imshow(plot_matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    # Annotate cells
    for i in range(len(all_types)):
        for j in range(len(vaf_bins)):
            v = matrix[i, j]
            if not math.isnan(v):
                text_color = "black" if 0.3 < v < 0.8 else "white"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=10, color=text_color, fontweight="bold")
            else:
                ax.text(j, i, "N/A", ha="center", va="center",
                        fontsize=8, color="grey")

    ax.set_xticks(range(len(bin_labels)))
    ax.set_xticklabels(bin_labels, fontsize=10)
    ax.set_yticks(range(len(all_types)))
    ax.set_yticklabels(all_types, fontsize=10)
    ax.set_xlabel("True VAF bin", fontsize=11)
    ax.set_ylabel("Variant type", fontsize=11)
    ax.set_title("Sensitivity Heatmap: Variant Type × VAF Range", fontsize=12)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Sensitivity", fontsize=10)

    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# 8. Threshold sweep (sensitivity/precision vs rVAF threshold)
# ---------------------------------------------------------------------------

def plot_threshold_sweep(
    sweep_df: pd.DataFrame,
    out_path: pathlib.Path,
) -> None:
    """Sensitivity and precision vs rVAF threshold. Two y-axes or two panels."""
    if sweep_df.empty:
        print("  Skipping: no threshold sweep data.")
        return

    sweep = sweep_df.dropna(subset=["vaf_threshold", "sensitivity"]).sort_values("vaf_threshold")
    if sweep.empty:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_SIZE_WIDE)

    # Panel 1: Sensitivity and precision vs threshold
    ax1.plot(sweep["vaf_threshold"], sweep["sensitivity"], color=KMERDET_COLOR,
             lw=2, label="Sensitivity")
    if "precision" in sweep.columns:
        ax1.plot(sweep["vaf_threshold"], sweep["precision"], color=THESIS_COLOR,
                 lw=2, ls="--", label="Precision")
    if "f1" in sweep.columns:
        f1_vals = sweep["f1"].apply(_safe_float)
        ax1.plot(sweep["vaf_threshold"], f1_vals, color="purple",
                 lw=1.5, ls=":", label="F1")

    ax1.set_xscale("log")
    ax1.set_xlabel("Minimum rVAF threshold", fontsize=10)
    ax1.set_ylabel("Score", fontsize=10)
    ax1.set_ylim(0, 1.05)
    ax1.set_title("Sensitivity / Precision vs rVAF Threshold", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)

    # Annotate clinical operating points
    for thresh, label in [(0.001, "0.1%"), (0.005, "0.5%"), (0.01, "1%")]:
        ax1.axvline(thresh, ls=":", color=NEUTRAL_COLOR, lw=0.8, alpha=0.6)
        ax1.text(thresh, 0.02, label, rotation=90, va="bottom", ha="right",
                 fontsize=7, color=NEUTRAL_COLOR)

    # Panel 2: ROC-like curve (sensitivity vs 1-specificity)
    if "specificity" in sweep.columns:
        fpr = 1 - sweep["specificity"].apply(lambda x: _safe_float(x) or 0)
        tpr = sweep["sensitivity"].apply(lambda x: _safe_float(x) or 0)
        ax2.plot(fpr, tpr, color=KMERDET_COLOR, lw=2, label="kmerdet")
        ax2.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Random")
        ax2.set_xlabel("1 - Specificity (FPR)", fontsize=10)
        ax2.set_ylabel("Sensitivity (TPR)", fontsize=10)
        ax2.set_title("ROC-like Curve (varying rVAF threshold)", fontsize=11)
        ax2.legend(fontsize=9)
        ax2.grid(alpha=0.3)
        ax2.set_xlim(-0.02, 1.02)
        ax2.set_ylim(-0.02, 1.05)
    else:
        ax2.text(0.5, 0.5, "No specificity data\navailable",
                 ha="center", va="center", fontsize=11, transform=ax2.transAxes)
        ax2.axis("off")

    fig.suptitle("rVAF Threshold Sweep Analysis", fontsize=12)
    fig.tight_layout()
    savefig(fig, out_path)


# ---------------------------------------------------------------------------
# Dispatcher: load all available data and produce all plots
# ---------------------------------------------------------------------------

def load_all_data(
    results_dir: Optional[pathlib.Path],
    metrics_json: Optional[pathlib.Path],
    per_type_tsv: Optional[pathlib.Path],
    per_vaf_tsv: Optional[pathlib.Path],
    per_variant_tsv: Optional[pathlib.Path],
    sweep_tsv: Optional[pathlib.Path],
) -> dict:
    """Load benchmark data from multiple possible sources into a single dict."""
    data: dict = {
        "metrics": {},
        "per_type": pd.DataFrame(),
        "per_vaf": pd.DataFrame(),
        "per_variant": pd.DataFrame(),
        "sweep": pd.DataFrame(),
        "timing": [],
    }

    # Load from metrics JSON
    if metrics_json and metrics_json.exists():
        with open(metrics_json) as f:
            data["metrics"] = json.load(f)

        # Extract embedded dataframes if present
        if "aggregate" in data["metrics"]:
            pass  # aggregated data, no per-type/per-vaf in top-level

    # Load from results_dir (auto-discover)
    if results_dir and results_dir.exists():
        if not metrics_json:
            mj = results_dir / "metrics.json"
            if mj.exists():
                with open(mj) as f:
                    data["metrics"] = json.load(f)

        # Performance timing records from metrics.json
        perf = data["metrics"].get("performance", {})
        if perf and "summary" in perf:
            # Convert to list of records for plotting
            data["timing"] = [
                {"label": k, **v}
                for k, v in perf.items()
                if isinstance(v, dict) and "wall_sec" in v
            ]

        # TSV files in results_dir
        for attr, filename in [
            ("per_type",    "per_type.tsv"),
            ("per_vaf",     "per_vaf_bin.tsv"),
            ("per_variant", "per_variant.tsv"),
            ("sweep",       "threshold_sweep.tsv"),
        ]:
            candidate = results_dir / filename
            if candidate.exists() and data[attr].empty:
                try:
                    data[attr] = pd.read_csv(candidate, sep="\t")
                except Exception as e:
                    print(f"  Warning: could not load {candidate}: {e}", file=sys.stderr)

    # Override with explicitly provided files
    for attr, path in [
        ("per_type",    per_type_tsv),
        ("per_vaf",     per_vaf_tsv),
        ("per_variant", per_variant_tsv),
        ("sweep",       sweep_tsv),
    ]:
        if path and path.exists():
            try:
                data[attr] = pd.read_csv(path, sep="\t")
            except Exception as e:
                print(f"  Warning: could not load {path}: {e}", file=sys.stderr)

    # Also check metrics.json for embedded dataframes
    if data["per_type"].empty and "per_type" in data["metrics"]:
        data["per_type"] = pd.DataFrame(data["metrics"]["per_type"])
    if data["per_vaf"].empty and "per_vaf_bin" in data["metrics"]:
        data["per_vaf"] = pd.DataFrame(data["metrics"]["per_vaf_bin"])
    if data["sweep"].empty and "threshold_sweep" in data["metrics"]:
        data["sweep"] = pd.DataFrame(data["metrics"]["threshold_sweep"])

    return data


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--results-dir",   type=pathlib.Path, help="analyze_results.py output directory")
    p.add_argument("--metrics",       type=pathlib.Path, help="metrics.json from analyze_results.py")
    p.add_argument("--per-type",      type=pathlib.Path, help="per_type.tsv")
    p.add_argument("--per-vaf",       type=pathlib.Path, help="per_vaf_bin.tsv")
    p.add_argument("--per-variant",   type=pathlib.Path, help="per_variant.tsv")
    p.add_argument("--sweep",         type=pathlib.Path, help="threshold_sweep.tsv")
    p.add_argument("--multi-k-per-type", type=pathlib.Path,
                   help="per_type.tsv from a multi-k run for comparison")
    p.add_argument("--output-dir",    type=pathlib.Path, default=pathlib.Path("figures"),
                   help="Directory to write figures [default: figures/]")
    p.add_argument("--format",        choices=["png", "pdf", "svg"], default="png",
                   help="Figure format [default: png]")
    p.add_argument("--dpi",           type=int, default=150, help="Figure DPI [default: 150]")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    global FIG_DPI
    FIG_DPI = args.dpi

    ext = args.format

    print(f"Loading benchmark data...")
    data = load_all_data(
        results_dir=args.results_dir,
        metrics_json=args.metrics,
        per_type_tsv=args.per_type,
        per_vaf_tsv=args.per_vaf,
        per_variant_tsv=args.per_variant,
        sweep_tsv=args.sweep,
    )

    print(f"\nGenerating figures -> {out_dir}/")

    # 1. Sensitivity by VAF
    plot_sensitivity_by_vaf(
        data["per_vaf"],
        out_dir / f"sensitivity_by_vaf.{ext}",
    )

    # 2. Runtime scaling
    plot_runtime_scaling(
        data["timing"],
        out_dir / f"runtime_scaling.{ext}",
    )

    # 3. Memory usage
    plot_memory_usage(
        data["timing"],
        out_dir / f"memory_by_dataset.{ext}",
    )

    # 4. Variant type breakdown
    plot_variant_type_breakdown(
        data["per_type"],
        out_dir / f"variant_type_breakdown.{ext}",
    )

    # 5. Multi-k comparison (only if second per-type file provided)
    if args.multi_k_per_type and args.multi_k_per_type.exists():
        multi_k_df = pd.read_csv(args.multi_k_per_type, sep="\t")
        plot_multi_k_comparison(
            data["per_type"],
            multi_k_df,
            out_dir / f"multi_k_comparison.{ext}",
        )
    else:
        print("  Skipping multi_k_comparison (no --multi-k-per-type provided).")

    # 6. Filtering Venn (need pre/post filter data)
    metrics = data["metrics"]
    confusion = metrics.get("confusion", {})
    if confusion.get("tp") is not None:
        # Rough estimate: pre-filter = tp + fp, post-filter = tp
        tp   = confusion.get("tp", 0)
        fp   = confusion.get("fp", 0)
        fn   = confusion.get("fn", confusion.get("fn_count", 0))
        n_before = tp + fp + fn  # estimated unfiltered
        n_after  = tp
        plot_filtering_venn(
            n_before=n_before,
            n_after=n_after,
            n_tp=tp,
            n_fp_removed=fp,
            n_fn_added=fn,
            out_path=out_dir / f"filtering_venn.{ext}",
        )
    else:
        print("  Skipping filtering_venn (no confusion matrix data).")

    # 7. Sensitivity heatmap
    plot_sensitivity_heatmap(
        data["per_variant"],
        out_dir / f"sensitivity_heatmap.{ext}",
    )

    # 8. Threshold sweep
    plot_threshold_sweep(
        data["sweep"],
        out_dir / f"threshold_sweep.{ext}",
    )

    print(f"\nAll figures saved to: {out_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
