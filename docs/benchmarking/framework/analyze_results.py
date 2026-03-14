#!/usr/bin/env python3
"""analyze_results.py — Compute performance metrics from kmerdet benchmark output.

Parses kmerdet TSV detection output and ground truth files, then computes:
  - Sensitivity, specificity, PPV (precision), NPV, F1, MCC
  - Per-variant-type breakdown
  - Per-VAF-bin breakdown (0-0.1%, 0.1-1%, 1-5%, 5-50%)
  - Confusion matrices
  - McNemar's test for kmerdet vs km comparison
  - ROC-like threshold sweep
  - Performance metrics (runtime, memory, throughput)

Usage:
    python3 analyze_results.py --results-dir /path/to/results
    python3 analyze_results.py --detected detected.tsv --truth ground_truth.tsv
    python3 analyze_results.py --benchmark-json benchmark.json

Output:
    metrics.json       — All computed metrics in JSON
    summary.tsv        — Human-readable summary table
    per_type.tsv       — Per-variant-type sensitivity breakdown
    per_vaf_bin.tsv    — Sensitivity by VAF bin
    per_variant.tsv    — Per-variant classification with rVAF vs true VAF
    comparison.tsv     — kmerdet vs km head-to-head (if both available)
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import confusion_matrix, matthews_corrcoef


# ---------------------------------------------------------------------------
# Constants: VAF bins matching thesis analysis
# ---------------------------------------------------------------------------

DEFAULT_VAF_BINS = [
    (0.0,    0.001,  "0–0.1%"),
    (0.001,  0.01,   "0.1–1%"),
    (0.01,   0.05,   "1–5%"),
    (0.05,   1.0,    "5–100%"),
]

ULTRA_FINE_VAF_BINS = [
    (0.0,       0.000001, "0–0.0001%"),
    (0.000001,  0.00001,  "0.0001–0.001%"),
    (0.00001,   0.0001,   "0.001–0.01%"),
    (0.0001,    0.001,    "0.01–0.1%"),
    (0.001,     0.01,     "0.1–1%"),
    (0.01,      0.05,     "1–5%"),
    (0.05,      1.0,      "5–100%"),
]

VARIANT_TYPE_ORDER = [
    "Substitution", "Insertion", "Deletion", "ITD", "Complex",
]

# Thesis baseline values for regression detection
THESIS_BASELINES = {
    "Substitution": 0.77,
    "Insertion":    0.38,
    "Deletion":     0.38,
    "ITD":          None,
    "overall_snv":  0.77,
    "overall_indel": 0.38,
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_detection_tsv(path: pathlib.Path) -> pd.DataFrame:
    """Load a kmerdet detection TSV file into a DataFrame.

    Expected columns (km-compatible format):
      sample, target, type, variant_name, rVAF, expression, min_coverage,
      start_kmer_count, ref_sequence, alt_sequence, info,
      chrom, pos, ref_allele, alt_allele, pvalue, qual, ci_lower, ci_upper
    """
    col_names = [
        "sample", "target", "variant_type", "variant_name",
        "rvaf", "expression", "min_coverage", "start_kmer_count",
        "ref_sequence", "alt_sequence", "info",
        "chrom", "pos", "ref_allele", "alt_allele",
        "pvalue", "qual", "ci_lower", "ci_upper",
    ]
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            names=col_names,
            header=0,      # skip the header row
            dtype=str,
            na_values=["", "NA", "N/A", "."],
        )
    except Exception as exc:
        raise ValueError(f"Failed to load detection TSV {path}: {exc}") from exc

    # Coerce numeric columns
    for col in ["rvaf", "expression", "min_coverage", "start_kmer_count",
                "pos", "pvalue", "qual", "ci_lower", "ci_upper"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def load_ground_truth(path: pathlib.Path) -> pd.DataFrame:
    """Load a ground truth TSV file.

    Expected columns:
      chrom, pos, ref, alt, type, true_vaf[, category]
    """
    col_names = ["chrom", "pos", "ref_allele", "alt_allele",
                 "variant_type", "true_vaf", "category"]
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            names=col_names,
            header=0,
            dtype=str,
            na_values=["", "NA", "N/A", "."],
        )
    except Exception as exc:
        raise ValueError(f"Failed to load ground truth {path}: {exc}") from exc

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["true_vaf"] = pd.to_numeric(df["true_vaf"], errors="coerce").fillna(0.0)

    # Normalize chromosome names (strip "chr" prefix for comparison)
    df["chrom_norm"] = df["chrom"].str.replace(r"^chr", "", regex=True).str.lower()

    return df


def load_benchmark_json(path: pathlib.Path) -> dict:
    """Load a pre-computed kmerdet benchmark JSON file."""
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Variant matching
# ---------------------------------------------------------------------------

def normalize_chrom(c: Optional[str]) -> Optional[str]:
    if c is None or (isinstance(c, float) and math.isnan(c)):
        return None
    return str(c).lower().removeprefix("chr")


def match_calls_to_truth(
    calls: pd.DataFrame,
    truth: pd.DataFrame,
) -> pd.DataFrame:
    """Match each detected call to a ground truth variant.

    Returns a copy of `truth` with additional columns:
      - detected (bool)
      - measured_rvaf (float | NaN)
      - vaf_error (float | NaN)
      - status ("TP" | "FN" | "FP_truth" | "TN")

    Also returns a list of FP calls (detected but not in truth).
    """
    truth = truth.copy()
    truth["detected"] = False
    truth["measured_rvaf"] = float("nan")
    truth["vaf_error"] = float("nan")
    truth["status"] = "FN"

    # Track which calls were matched
    calls = calls.copy()
    calls["matched"] = False

    # Normalize chromosomes in calls
    calls["chrom_norm"] = calls["chrom"].apply(normalize_chrom)

    for idx, tv in truth.iterrows():
        chrom_t = normalize_chrom(tv["chrom"])
        pos_t = tv["pos"]
        ref_t = str(tv["ref_allele"]) if pd.notna(tv.get("ref_allele")) else None
        alt_t = str(tv["alt_allele"]) if pd.notna(tv.get("alt_allele")) else None

        # Find matching calls
        mask = (
            (calls["chrom_norm"] == chrom_t) &
            (calls["pos"] == pos_t) &
            (calls["ref_allele"].astype(str) == ref_t) &
            (calls["alt_allele"].astype(str) == alt_t) &
            (~calls["matched"])
        )
        matching = calls[mask]

        if not matching.empty:
            # Prefer highest rVAF if multiple matches
            best_idx = matching["rvaf"].idxmax()
            calls.at[best_idx, "matched"] = True
            measured = calls.at[best_idx, "rvaf"]
            truth.at[idx, "detected"] = True
            truth.at[idx, "measured_rvaf"] = measured
            if tv["true_vaf"] > 0:
                truth.at[idx, "vaf_error"] = measured - tv["true_vaf"]
                truth.at[idx, "status"] = "TP"
            else:
                truth.at[idx, "status"] = "FP_truth"
        else:
            if tv["true_vaf"] > 0:
                truth.at[idx, "status"] = "FN"
            else:
                truth.at[idx, "status"] = "TN"

    # Unmatched non-Reference calls are FPs
    ref_types = {"Reference", "reference"}
    fps = calls[
        (~calls["matched"]) &
        (~calls["variant_type"].isin(ref_types))
    ].copy()
    fps["status"] = "FP_call"
    fps["true_vaf"] = 0.0

    return truth, fps


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def compute_confusion(truth_matched: pd.DataFrame, fps: pd.DataFrame) -> dict:
    """Compute TP, FP, FN, TN from matched truth and unmatched FP calls."""
    tp = int((truth_matched["status"] == "TP").sum())
    fn = int((truth_matched["status"] == "FN").sum())
    tn = int((truth_matched["status"] == "TN").sum())
    fp_truth = int((truth_matched["status"] == "FP_truth").sum())
    fp_calls = len(fps)
    fp = fp_truth + fp_calls
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "fp_truth": fp_truth, "fp_calls": fp_calls}


def safe_div(numerator: int, denominator: int) -> float:
    return float("nan") if denominator == 0 else numerator / denominator


def compute_metrics(cm: dict) -> dict:
    """Compute all standard binary classification metrics from a confusion dict."""
    tp, fp, fn, tn = cm["tp"], cm["fp"], cm["fn"], cm["tn"]
    sensitivity = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    precision   = safe_div(tp, tp + fp)
    npv         = safe_div(tn, tn + fn)
    f1          = (
        float("nan")
        if math.isnan(precision) or math.isnan(sensitivity) or (precision + sensitivity) == 0
        else 2 * precision * sensitivity / (precision + sensitivity)
    )
    accuracy    = safe_div(tp + tn, tp + tn + fp + fn)

    # Matthews Correlation Coefficient
    denom_sq = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = float("nan") if denom_sq == 0 else (tp * tn - fp * fn) / denom_sq

    return {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision":   precision,
        "npv":         npv,
        "f1":          f1,
        "accuracy":    accuracy,
        "mcc":         mcc,
        "tp":          tp, "fp": fp, "fn": fn, "tn": tn,
    }


def compute_per_type(truth_matched: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity breakdown by variant type (present truth variants only)."""
    rows = []
    present = truth_matched[truth_matched["true_vaf"] > 0].copy()

    for vtype in VARIANT_TYPE_ORDER + list(
        set(present["variant_type"].dropna().unique()) - set(VARIANT_TYPE_ORDER)
    ):
        sub = present[present["variant_type"] == vtype]
        if sub.empty:
            continue
        n_present = len(sub)
        n_tp = int((sub["status"] == "TP").sum())
        n_fn = int((sub["status"] == "FN").sum())
        sensitivity = safe_div(n_tp, n_tp + n_fn)
        baseline = THESIS_BASELINES.get(vtype)
        rows.append({
            "variant_type": vtype,
            "present":      n_present,
            "tp":           n_tp,
            "fn":           n_fn,
            "sensitivity":  sensitivity,
            "thesis_baseline": baseline,
            "vs_baseline": (
                sensitivity - baseline if baseline is not None and not math.isnan(sensitivity)
                else float("nan")
            ),
        })

    return pd.DataFrame(rows)


def compute_per_vaf_bin(
    truth_matched: pd.DataFrame,
    vaf_bins: list[tuple[float, float, str]] = DEFAULT_VAF_BINS,
) -> pd.DataFrame:
    """Sensitivity breakdown by true VAF bin."""
    present = truth_matched[truth_matched["true_vaf"] > 0]
    rows = []

    for low, high, label in vaf_bins:
        in_bin = present[
            (present["true_vaf"] >= low) &
            (present["true_vaf"] < high if high < 1.0 else present["true_vaf"] <= high)
        ]
        n = len(in_bin)
        n_tp = int((in_bin["status"] == "TP").sum())
        n_fn = int((in_bin["status"] == "FN").sum())
        rows.append({
            "bin_label":   label,
            "vaf_low":     low,
            "vaf_high":    high,
            "present":     n,
            "tp":          n_tp,
            "fn":          n_fn,
            "sensitivity": safe_div(n_tp, n_tp + n_fn),
        })

    return pd.DataFrame(rows)


def compute_per_variant(truth_matched: pd.DataFrame, fps: pd.DataFrame) -> pd.DataFrame:
    """Per-variant table with true VAF, measured rVAF, and status."""
    keep_cols = [
        "chrom", "pos", "ref_allele", "alt_allele", "variant_type",
        "true_vaf", "measured_rvaf", "vaf_error", "status", "category",
    ]
    existing = [c for c in keep_cols if c in truth_matched.columns]
    result = truth_matched[existing].copy()

    if not fps.empty:
        fp_rows = fps[["chrom", "pos", "ref_allele", "alt_allele",
                        "variant_type", "rvaf"]].copy()
        fp_rows = fp_rows.rename(columns={"rvaf": "measured_rvaf"})
        fp_rows["true_vaf"] = 0.0
        fp_rows["vaf_error"] = float("nan")
        fp_rows["status"] = "FP_call"
        result = pd.concat([result, fp_rows], ignore_index=True)

    return result


def compute_lod(
    per_vaf_df: pd.DataFrame,
    target_sensitivities: list[float] = [0.5, 0.8],
) -> dict:
    """Interpolate the sensitivity curve to find LOD (Limit of Detection).

    For each target sensitivity (e.g. 50%, 80%), finds the VAF at which the
    sensitivity curve crosses that threshold by linear interpolation between
    VAF bin midpoints.

    Args:
        per_vaf_df: DataFrame from compute_per_vaf_bin with columns
            vaf_low, vaf_high, sensitivity, present.
        target_sensitivities: List of sensitivity thresholds to find LOD for.

    Returns:
        Dict with keys like 'lod50', 'lod80' mapping to the interpolated VAF,
        or None if the threshold is never crossed.
    """
    if per_vaf_df.empty:
        return {f"lod{int(t*100)}": None for t in target_sensitivities}

    # Filter bins with actual data and compute midpoints
    valid = per_vaf_df[per_vaf_df["present"] > 0].copy()
    if valid.empty:
        return {f"lod{int(t*100)}": None for t in target_sensitivities}

    valid = valid.sort_values("vaf_low")
    midpoints = ((valid["vaf_low"] + valid["vaf_high"]) / 2).values
    sensitivities = valid["sensitivity"].values

    result = {}
    for target in target_sensitivities:
        key = f"lod{int(target * 100)}"
        lod_val = None

        # Walk from low VAF to high VAF; find first crossing of the target
        for i in range(len(sensitivities) - 1):
            s_low = sensitivities[i]
            s_high = sensitivities[i + 1]

            # Check if this segment crosses the target (low->high or high->low)
            if (math.isnan(s_low) or math.isnan(s_high)):
                continue

            # We want the VAF where sensitivity first reaches the target
            # going from low VAF (harder) to high VAF (easier)
            if s_low < target <= s_high:
                # Linear interpolation
                frac = (target - s_low) / (s_high - s_low) if s_high != s_low else 0.0
                lod_val = midpoints[i] + frac * (midpoints[i + 1] - midpoints[i])
                break
            elif s_low >= target:
                # Already above target at this bin — LOD is at or below this midpoint
                lod_val = midpoints[i]
                break

        result[key] = lod_val

    return result


def compute_threshold_sweep(
    calls: pd.DataFrame,
    truth: pd.DataFrame,
    thresholds: Optional[list[float]] = None,
) -> pd.DataFrame:
    """Compute (sensitivity, specificity, precision) at each rVAF threshold."""
    if thresholds is None:
        thresholds = np.logspace(-4, 0, 50).tolist() + [0.0]
        thresholds = sorted(set(thresholds))

    rows = []
    for thresh in thresholds:
        filtered = calls[
            (calls["rvaf"] >= thresh) &
            (~calls["variant_type"].isin({"Reference", "reference"}))
        ]
        truth_m, fps_m = match_calls_to_truth(filtered, truth)
        cm = compute_confusion(truth_m, fps_m)
        m = compute_metrics(cm)
        rows.append({
            "vaf_threshold": thresh,
            "sensitivity":   m["sensitivity"],
            "specificity":   m["specificity"],
            "precision":     m["precision"],
            "f1":            m["f1"],
            "tp":            m["tp"],
            "fp":            m["fp"],
            "fn":            m["fn"],
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# McNemar's test for kmerdet vs km comparison
# ---------------------------------------------------------------------------

def mcnemar_test(
    kmerdet_tp_mask: np.ndarray,
    km_tp_mask: np.ndarray,
) -> dict:
    """McNemar's test comparing two detectors on the same set of truth variants.

    Both arrays must be boolean with length = number of truth positives.
    Returns the test statistic, p-value, and a qualitative interpretation.
    """
    # Contingency table:
    # b = kmerdet correct, km incorrect
    # c = kmerdet incorrect, km correct
    b = int(( kmerdet_tp_mask & ~km_tp_mask).sum())
    c = int((~kmerdet_tp_mask &  km_tp_mask).sum())

    # McNemar's test (with continuity correction for small b+c)
    n = b + c
    if n == 0:
        return {"b": 0, "c": 0, "statistic": 0.0, "pvalue": 1.0,
                "interpretation": "No discordant pairs; detectors are identical on this data."}

    # Chi-square with continuity correction
    chi2 = (abs(b - c) - 1) ** 2 / n if n >= 20 else (abs(b - c)) ** 2 / n
    pvalue = float(stats.chi2.sf(chi2, df=1))

    if pvalue < 0.001:
        interp = f"Significant difference (p={pvalue:.4f}); "
        interp += "kmerdet is better." if b > c else "km is better."
    elif pvalue < 0.05:
        interp = f"Marginally significant (p={pvalue:.4f})."
    else:
        interp = f"No significant difference (p={pvalue:.4f})."

    return {
        "b_kmerdet_only": b,
        "c_km_only": c,
        "statistic": chi2,
        "pvalue": pvalue,
        "interpretation": interp,
    }


# ---------------------------------------------------------------------------
# Performance analysis
# ---------------------------------------------------------------------------

def parse_timing_files(results_dir: pathlib.Path) -> pd.DataFrame:
    """Parse all *_timing.json files in results_dir/performance/."""
    perf_dir = results_dir / "performance"
    if not perf_dir.exists():
        return pd.DataFrame()

    rows = []
    for f in sorted(perf_dir.glob("*_timing*.json")):
        try:
            with open(f) as fh:
                d = json.load(fh)
            rows.append({
                "file":        f.stem,
                "label":       d.get("label", f.stem),
                "wall_sec":    d.get("wall_sec"),
                "peak_rss_kb": d.get("peak_rss_kb"),
                "peak_rss_mb": (d["peak_rss_kb"] / 1024
                                if d.get("peak_rss_kb") is not None else None),
            })
        except Exception as e:
            print(f"Warning: could not parse {f}: {e}", file=sys.stderr)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def analyze_performance(timing_df: pd.DataFrame) -> dict:
    """Summarize performance metrics from timing data."""
    if timing_df.empty:
        return {"error": "No timing data available."}

    result = {}
    for _, row in timing_df.iterrows():
        result[row["label"]] = {
            "wall_sec":    row["wall_sec"],
            "peak_rss_mb": row.get("peak_rss_mb"),
        }

    # Find the main detection timing
    detect_rows = timing_df[timing_df["label"].str.startswith("detect_")]
    if not detect_rows.empty:
        main = detect_rows.iloc[0]
        result["summary"] = {
            "fastest_detect_sec": detect_rows["wall_sec"].min(),
            "slowest_detect_sec": detect_rows["wall_sec"].max(),
            "vs_thesis_6min": (
                detect_rows["wall_sec"].min() / 360.0
                if detect_rows["wall_sec"].min() is not None else None
            ),
        }

    return result


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def print_summary_table(
    metrics: dict,
    per_type: pd.DataFrame,
    per_vaf_bin: pd.DataFrame,
) -> None:
    """Print a human-readable summary to stdout."""
    print("\n=== kmerdet Benchmark Results ===\n")

    # Overall metrics
    cm = metrics.get("confusion", {})
    m  = metrics.get("overall", {})
    print(f"Overall:")
    print(f"  Sensitivity:  {m.get('sensitivity', float('nan')):.3f}")
    print(f"  Specificity:  {m.get('specificity', float('nan')):.3f}")
    print(f"  Precision:    {m.get('precision',   float('nan')):.3f}")
    print(f"  NPV:          {m.get('npv',         float('nan')):.3f}")
    print(f"  F1:           {m.get('f1',           float('nan')):.3f}")
    print(f"  MCC:          {m.get('mcc',          float('nan')):.3f}")
    print(f"  Accuracy:     {m.get('accuracy',     float('nan')):.3f}")
    print(f"  TP={cm.get('tp',0)}  FP={cm.get('fp',0)}  "
          f"FN={cm.get('fn',0)}  TN={cm.get('tn',0)}")
    print()

    # Per-type table
    if not per_type.empty:
        print("Per-type sensitivity:")
        print(f"  {'Type':<20} {'Present':>8} {'TP':>5} {'FN':>5} "
              f"{'Sensitivity':>12} {'vs Thesis':>10}")
        print("  " + "-" * 64)
        for _, row in per_type.iterrows():
            sens = f"{row['sensitivity']:.3f}" if not math.isnan(row['sensitivity']) else " N/A"
            bl   = (f"{row['vs_baseline']:+.3f}"
                    if not math.isnan(row.get('vs_baseline', float('nan'))) else "   N/A")
            print(f"  {row['variant_type']:<20} {int(row['present']):>8} "
                  f"{int(row['tp']):>5} {int(row['fn']):>5} {sens:>12} {bl:>10}")
    print()

    # Per-VAF-bin table
    if not per_vaf_bin.empty:
        print("Sensitivity by VAF bin:")
        print(f"  {'VAF range':<15} {'Present':>8} {'TP':>5} {'FN':>5} {'Sensitivity':>12}")
        print("  " + "-" * 44)
        for _, row in per_vaf_bin.iterrows():
            if row["present"] == 0:
                continue
            sens = f"{row['sensitivity']:.3f}" if not math.isnan(row['sensitivity']) else " N/A"
            print(f"  {row['bin_label']:<15} {int(row['present']):>8} "
                  f"{int(row['tp']):>5} {int(row['fn']):>5} {sens:>12}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = p.add_mutually_exclusive_group()
    source.add_argument(
        "--results-dir",
        type=pathlib.Path,
        help="Directory produced by run_benchmarks.sh. "
             "Auto-discovers all accuracy/*.json and performance/ files.",
    )
    source.add_argument(
        "--detected",
        type=pathlib.Path,
        help="kmerdet detection TSV output file.",
    )
    source.add_argument(
        "--benchmark-json",
        type=pathlib.Path,
        help="Pre-computed kmerdet benchmark JSON (from kmerdet benchmark --format json).",
    )

    p.add_argument(
        "--truth",
        type=pathlib.Path,
        help="Ground truth TSV file (required with --detected).",
    )
    p.add_argument(
        "--km-detected",
        type=pathlib.Path,
        help="km detection TSV for head-to-head comparison with McNemar's test.",
    )
    p.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path("analysis_output"),
        help="Directory to write output files [default: analysis_output/].",
    )
    p.add_argument(
        "--sweep",
        action="store_true",
        help="Run VAF threshold sweep (generates ROC-like curve data).",
    )
    p.add_argument(
        "--vaf-bins",
        type=str,
        default=None,
        help="VAF bin boundaries as comma-separated floats (e.g. '0,0.001,0.01,0.1,1.0'). "
             "Each adjacent pair defines one bin.",
    )
    p.add_argument(
        "--no-summary",
        action="store_true",
        help="Suppress printed summary table.",
    )
    p.add_argument(
        "--ultra-fine",
        action="store_true",
        help="Use ultra-fine VAF bins for sensitivity-first analysis "
             "(extends into sub-0.1%% regime).",
    )

    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------- Load data -------
    all_metrics: dict = {}
    per_type_df = pd.DataFrame()
    per_vaf_df  = pd.DataFrame()
    per_variant_df = pd.DataFrame()
    timing_df = pd.DataFrame()

    if args.vaf_bins:
        boundaries = [float(x) for x in args.vaf_bins.split(",")]
        custom_bins = [
            (boundaries[i], boundaries[i+1],
             f"{boundaries[i]*100:.3g}–{boundaries[i+1]*100:.3g}%")
            for i in range(len(boundaries) - 1)
        ]
    elif args.ultra_fine:
        custom_bins = ULTRA_FINE_VAF_BINS
    else:
        custom_bins = DEFAULT_VAF_BINS

    if args.benchmark_json:
        # Already-computed benchmark JSON from kmerdet benchmark subcommand
        print(f"Loading pre-computed benchmark: {args.benchmark_json}")
        raw = load_benchmark_json(args.benchmark_json)
        cm  = raw.get("confusion", {})
        s   = raw.get("summary", {})
        all_metrics["confusion"] = cm
        all_metrics["overall"] = {
            "sensitivity": s.get("sensitivity"),
            "specificity": s.get("specificity"),
            "precision":   s.get("precision"),
            "npv":         s.get("npv"),
            "f1":          s.get("f1"),
            "accuracy":    s.get("accuracy"),
            "mcc":         float("nan"),  # not in benchmark JSON
        }
        per_type_df = pd.DataFrame(raw.get("per_type", []))
        per_vaf_df  = pd.DataFrame(raw.get("per_vaf_bin", []))
        per_variant_df = pd.DataFrame(raw.get("per_variant", []))

    elif args.detected:
        if args.truth is None:
            print("ERROR: --truth is required when using --detected", file=sys.stderr)
            return 1

        print(f"Loading detection results: {args.detected}")
        print(f"Loading ground truth:      {args.truth}")

        calls = load_detection_tsv(args.detected)
        truth = load_ground_truth(args.truth)

        print(f"  Calls: {len(calls)} rows  (non-ref: "
              f"{(~calls['variant_type'].isin({'Reference','reference'})).sum()})")
        print(f"  Truth: {len(truth)} variants  "
              f"(present: {(truth['true_vaf'] > 0).sum()})")

        truth_matched, fps = match_calls_to_truth(calls, truth)
        cm = compute_confusion(truth_matched, fps)
        overall_metrics = compute_metrics(cm)
        all_metrics["confusion"] = cm
        all_metrics["overall"]   = overall_metrics
        per_type_df    = compute_per_type(truth_matched)
        per_vaf_df     = compute_per_vaf_bin(truth_matched, custom_bins)
        per_variant_df = compute_per_variant(truth_matched, fps)

        if args.sweep:
            print("Running threshold sweep...")
            sweep_df = compute_threshold_sweep(calls, truth)
            sweep_df.to_csv(out_dir / "threshold_sweep.tsv", sep="\t", index=False)
            all_metrics["threshold_sweep"] = sweep_df.to_dict(orient="records")

        # McNemar's test if km results are provided
        if args.km_detected:
            print(f"Loading km comparison results: {args.km_detected}")
            km_calls = load_detection_tsv(args.km_detected)
            km_truth_m, km_fps = match_calls_to_truth(km_calls, truth)

            present = truth_matched[truth_matched["true_vaf"] > 0]
            km_present = km_truth_m[km_truth_m["true_vaf"] > 0]

            kmerdet_tp = (present["status"] == "TP").to_numpy()
            km_tp = (km_present["status"] == "TP").to_numpy()

            if len(kmerdet_tp) == len(km_tp):
                mcnemar = mcnemar_test(kmerdet_tp, km_tp)
                all_metrics["mcnemar"] = mcnemar
                print(f"\nMcNemar's test (kmerdet vs km):")
                print(f"  b (kmerdet only): {mcnemar['b_kmerdet_only']}")
                print(f"  c (km only):      {mcnemar['c_km_only']}")
                print(f"  p-value:          {mcnemar['pvalue']:.4f}")
                print(f"  {mcnemar['interpretation']}")

    elif args.results_dir:
        # Aggregate all benchmark JSON files from run_benchmarks.sh output
        print(f"Loading results from directory: {args.results_dir}")
        acc_dir = args.results_dir / "accuracy"
        if not acc_dir.exists():
            print(f"ERROR: No accuracy/ subdirectory in {args.results_dir}", file=sys.stderr)
            return 1

        benchmark_files = sorted(acc_dir.glob("*_benchmark.json"))
        if not benchmark_files:
            print(f"ERROR: No *_benchmark.json files in {acc_dir}", file=sys.stderr)
            return 1

        aggregate_rows = []
        for bf in benchmark_files:
            raw = load_benchmark_json(bf)
            cm  = raw.get("confusion", {})
            s   = raw.get("summary", {})
            aggregate_rows.append({
                "dataset":     bf.stem.replace("_benchmark", ""),
                "sensitivity": s.get("sensitivity"),
                "specificity": s.get("specificity"),
                "precision":   s.get("precision"),
                "f1":          s.get("f1"),
                "tp":          cm.get("tp", 0),
                "fp":          cm.get("fp", 0),
                "fn":          cm.get("fn_count", 0),
                "tn":          cm.get("tn", 0),
            })

        agg_df = pd.DataFrame(aggregate_rows)
        agg_df.to_csv(out_dir / "aggregate_accuracy.tsv", sep="\t", index=False)
        all_metrics["aggregate"] = aggregate_rows
        print(f"Aggregated {len(aggregate_rows)} benchmark runs.")

        # Use the first benchmark as the "main" for detailed breakdown
        if benchmark_files:
            raw = load_benchmark_json(benchmark_files[0])
            all_metrics["confusion"] = raw.get("confusion", {})
            all_metrics["overall"]   = raw.get("summary", {})
            per_type_df    = pd.DataFrame(raw.get("per_type", []))
            per_vaf_df     = pd.DataFrame(raw.get("per_vaf_bin", []))
            per_variant_df = pd.DataFrame(raw.get("per_variant", []))

        timing_df = parse_timing_files(args.results_dir)

    else:
        print("ERROR: provide --results-dir, --detected, or --benchmark-json", file=sys.stderr)
        return 1

    # ------- LOD computation -------
    if not per_vaf_df.empty:
        lod_metrics = compute_lod(per_vaf_df)
        all_metrics["lod"] = lod_metrics
        lod_strs = []
        for key, val in lod_metrics.items():
            if val is not None:
                lod_strs.append(f"{key}={val:.6f}")
            else:
                lod_strs.append(f"{key}=N/A")
        print(f"LOD metrics: {', '.join(lod_strs)}")

    # ------- Performance analysis -------
    if not timing_df.empty:
        perf = analyze_performance(timing_df)
        all_metrics["performance"] = perf
        timing_df.to_csv(out_dir / "performance.tsv", sep="\t", index=False)

    # ------- Write outputs -------
    metrics_path = out_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(all_metrics, f, indent=2, default=lambda x: None if (
            isinstance(x, float) and math.isnan(x)) else x)
    print(f"\nMetrics written: {metrics_path}")

    summary_path = out_dir / "summary.tsv"
    with open(summary_path, "w") as f:
        m = all_metrics.get("overall", {})
        cm = all_metrics.get("confusion", {})
        f.write("metric\tvalue\n")
        for key in ["sensitivity", "specificity", "precision", "npv",
                    "f1", "mcc", "accuracy"]:
            v = m.get(key)
            f.write(f"{key}\t{v:.4f}\n" if v is not None and not (
                isinstance(v, float) and math.isnan(v)) else f"{key}\tNA\n")
        f.write(f"tp\t{cm.get('tp', 0)}\n")
        f.write(f"fp\t{cm.get('fp', 0)}\n")
        f.write(f"fn\t{cm.get('fn', cm.get('fn_count', 0))}\n")
        f.write(f"tn\t{cm.get('tn', 0)}\n")
    print(f"Summary written: {summary_path}")

    if not per_type_df.empty:
        per_type_df.to_csv(out_dir / "per_type.tsv", sep="\t", index=False)
        print(f"Per-type:  {out_dir / 'per_type.tsv'}")

    if not per_vaf_df.empty:
        per_vaf_df.to_csv(out_dir / "per_vaf_bin.tsv", sep="\t", index=False)
        print(f"Per-VAF:   {out_dir / 'per_vaf_bin.tsv'}")

    if not per_variant_df.empty:
        per_variant_df.to_csv(out_dir / "per_variant.tsv", sep="\t", index=False)
        print(f"Per-var:   {out_dir / 'per_variant.tsv'}")

    # ------- Print summary -------
    if not args.no_summary:
        print_summary_table(all_metrics, per_type_df, per_vaf_df)

    return 0


if __name__ == "__main__":
    sys.exit(main())
