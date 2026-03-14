#!/usr/bin/env python3
"""compare_results.py — Side-by-side analysis of km vs kmerdet benchmark results.

Loads benchmark JSON files from both tools, builds a paired concordance matrix,
computes per-tool sensitivity/specificity/F1, and runs McNemar's test for
statistical significance.

Outputs:
  comparison_summary.tsv   — Per-tool metrics side-by-side
  concordance_matrix.tsv   — Variant-level concordance (both TP, both FN, etc.)
  comparison_report.md     — Human-readable report with interpretations

Usage:
    python3 compare_results.py \\
        --km-dir results/km \\
        --kmerdet-dir results/kmerdet \\
        --truth ground_truth.tsv \\
        --output-dir results/comparison
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


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VARIANT_TYPE_ORDER = [
    "Substitution", "Insertion", "Deletion", "ITD", "Complex",
]

DEFAULT_VAF_BINS = [
    (0.0,    0.001,  "0-0.1%"),
    (0.001,  0.01,   "0.1-1%"),
    (0.01,   0.05,   "1-5%"),
    (0.05,   1.0,    "5-100%"),
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_benchmark_json(path: pathlib.Path) -> dict:
    with open(path) as f:
        return json.load(f)


def discover_benchmark_files(tool_dir: pathlib.Path) -> list[tuple[str, pathlib.Path]]:
    """Find all *_benchmark_*.json files and extract parameter suffix."""
    results = []
    for f in sorted(tool_dir.glob("*_benchmark_*.json")):
        # Extract suffix like "c2_r0.05" from filename
        stem = f.stem
        # Pattern: tool_benchmark_SUFFIX
        parts = stem.split("_benchmark_")
        suffix = parts[1] if len(parts) > 1 else stem
        results.append((suffix, f))
    return results


def safe_div(num: int, denom: int) -> float:
    return float("nan") if denom == 0 else num / denom


# ---------------------------------------------------------------------------
# Metrics extraction
# ---------------------------------------------------------------------------

def extract_metrics(bench: dict) -> dict:
    """Extract key metrics from a benchmark JSON."""
    s = bench.get("summary", {})
    c = bench.get("confusion", {})
    return {
        "sensitivity": s.get("sensitivity"),
        "specificity": s.get("specificity"),
        "precision":   s.get("precision"),
        "f1":          s.get("f1"),
        "tp":          c.get("tp", 0),
        "fp":          c.get("fp", 0),
        "fn":          c.get("fn_count", c.get("fn", 0)),
        "tn":          c.get("tn", 0),
    }


def extract_per_type(bench: dict) -> pd.DataFrame:
    """Extract per-type breakdown from benchmark JSON."""
    rows = bench.get("per_type", [])
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def extract_per_vaf(bench: dict) -> pd.DataFrame:
    """Extract per-VAF-bin breakdown from benchmark JSON."""
    rows = bench.get("per_vaf_bin", [])
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def extract_per_variant(bench: dict) -> pd.DataFrame:
    """Extract per-variant results from benchmark JSON."""
    rows = bench.get("per_variant", [])
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# McNemar's test
# ---------------------------------------------------------------------------

def mcnemar_test(
    kmerdet_tp: np.ndarray,
    km_tp: np.ndarray,
) -> dict:
    """McNemar's test comparing two detectors on the same truth set.

    Both arrays are boolean: True = variant detected (TP), False = missed (FN).
    """
    b = int(( kmerdet_tp & ~km_tp).sum())  # kmerdet only
    c = int((~kmerdet_tp &  km_tp).sum())  # km only

    n = b + c
    if n == 0:
        return {
            "b_kmerdet_only": 0, "c_km_only": 0,
            "statistic": 0.0, "pvalue": 1.0,
            "interpretation": "No discordant pairs; tools agree on all variants.",
        }

    # Chi-square with continuity correction for n >= 20
    if n >= 20:
        chi2 = (abs(b - c) - 1) ** 2 / n
    else:
        chi2 = (abs(b - c)) ** 2 / n
    pvalue = float(stats.chi2.sf(chi2, df=1))

    if pvalue < 0.001:
        interp = f"Highly significant (p={pvalue:.4e}); "
        interp += "kmerdet detects more." if b > c else "km detects more."
    elif pvalue < 0.05:
        interp = f"Significant (p={pvalue:.4f}); "
        interp += "kmerdet detects more." if b > c else "km detects more."
    else:
        interp = f"No significant difference (p={pvalue:.4f})."

    return {
        "b_kmerdet_only": b,
        "c_km_only": c,
        "both_detected": int((kmerdet_tp & km_tp).sum()),
        "both_missed": int((~kmerdet_tp & ~km_tp).sum()),
        "statistic": chi2,
        "pvalue": pvalue,
        "interpretation": interp,
    }


# ---------------------------------------------------------------------------
# Concordance matrix
# ---------------------------------------------------------------------------

def build_concordance(
    km_per_variant: pd.DataFrame,
    kmerdet_per_variant: pd.DataFrame,
) -> pd.DataFrame:
    """Build a variant-level concordance matrix between km and kmerdet.

    Matches on target name (or chrom+pos+alleles for real data).
    """
    if km_per_variant.empty or kmerdet_per_variant.empty:
        return pd.DataFrame()

    # For simulated data, match by variant identity
    # Use chrom + pos + ref_allele + alt_allele as key
    key_cols = ["chrom", "pos", "ref_allele", "alt_allele"]

    # Ensure key columns exist in both
    for col in key_cols:
        if col not in km_per_variant.columns or col not in kmerdet_per_variant.columns:
            return pd.DataFrame()

    km_v = km_per_variant.copy()
    kmerdet_v = kmerdet_per_variant.copy()

    # Create composite key
    km_v["key"] = km_v[key_cols].astype(str).agg("_".join, axis=1)
    kmerdet_v["key"] = kmerdet_v[key_cols].astype(str).agg("_".join, axis=1)

    # Merge
    merged = km_v[["key", "status", "variant_type", "true_vaf"]].rename(
        columns={"status": "km_status"}
    ).merge(
        kmerdet_v[["key", "status"]].rename(columns={"status": "kmerdet_status"}),
        on="key",
        how="outer",
    )

    # Classify concordance
    def concordance_class(row):
        km = row.get("km_status", "")
        kd = row.get("kmerdet_status", "")
        if km == "TP" and kd == "TP":
            return "both_TP"
        elif km == "FN" and kd == "FN":
            return "both_FN"
        elif km == "TP" and kd == "FN":
            return "km_only"
        elif km == "FN" and kd == "TP":
            return "kmerdet_only"
        elif km == "TN" and kd == "TN":
            return "both_TN"
        else:
            return "other"

    merged["concordance"] = merged.apply(concordance_class, axis=1)
    return merged


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_comparison_summary(
    km_metrics: dict,
    kmerdet_metrics: dict,
    km_per_type: pd.DataFrame,
    kmerdet_per_type: pd.DataFrame,
    km_per_vaf: pd.DataFrame,
    kmerdet_per_vaf: pd.DataFrame,
    param_suffix: str,
) -> pd.DataFrame:
    """Build a side-by-side summary DataFrame."""
    rows = []

    # Overall metrics
    for metric in ["sensitivity", "specificity", "precision", "f1"]:
        km_val = km_metrics.get(metric)
        kd_val = kmerdet_metrics.get(metric)
        diff = None
        if km_val is not None and kd_val is not None:
            try:
                diff = float(kd_val) - float(km_val)
            except (TypeError, ValueError):
                pass
        rows.append({
            "params": param_suffix,
            "category": "overall",
            "metric": metric,
            "km": km_val,
            "kmerdet": kd_val,
            "diff": diff,
        })

    # Confusion matrix
    for metric in ["tp", "fp", "fn", "tn"]:
        rows.append({
            "params": param_suffix,
            "category": "confusion",
            "metric": metric,
            "km": km_metrics.get(metric, 0),
            "kmerdet": kmerdet_metrics.get(metric, 0),
            "diff": None,
        })

    # Per-type sensitivity
    all_types = set()
    if not km_per_type.empty:
        all_types |= set(km_per_type["variant_type"])
    if not kmerdet_per_type.empty:
        all_types |= set(kmerdet_per_type["variant_type"])

    for vtype in VARIANT_TYPE_ORDER:
        if vtype not in all_types:
            continue
        km_sens = None
        kd_sens = None
        if not km_per_type.empty:
            row = km_per_type[km_per_type["variant_type"] == vtype]
            if not row.empty:
                km_sens = row.iloc[0].get("sensitivity")
        if not kmerdet_per_type.empty:
            row = kmerdet_per_type[kmerdet_per_type["variant_type"] == vtype]
            if not row.empty:
                kd_sens = row.iloc[0].get("sensitivity")
        diff = None
        if km_sens is not None and kd_sens is not None:
            try:
                diff = float(kd_sens) - float(km_sens)
            except (TypeError, ValueError):
                pass
        rows.append({
            "params": param_suffix,
            "category": f"type_{vtype}",
            "metric": "sensitivity",
            "km": km_sens,
            "kmerdet": kd_sens,
            "diff": diff,
        })

    # Per-VAF-bin sensitivity
    if not km_per_vaf.empty and not kmerdet_per_vaf.empty:
        for _, km_row in km_per_vaf.iterrows():
            label = km_row.get("bin_label", "")
            kd_row = kmerdet_per_vaf[kmerdet_per_vaf["bin_label"] == label]
            km_sens = km_row.get("sensitivity")
            kd_sens = kd_row.iloc[0]["sensitivity"] if not kd_row.empty else None
            diff = None
            if km_sens is not None and kd_sens is not None:
                try:
                    diff = float(kd_sens) - float(km_sens)
                except (TypeError, ValueError):
                    pass
            rows.append({
                "params": param_suffix,
                "category": f"vaf_{label}",
                "metric": "sensitivity",
                "km": km_sens,
                "kmerdet": kd_sens,
                "diff": diff,
            })

    return pd.DataFrame(rows)


def generate_report(
    summary_df: pd.DataFrame,
    concordance_df: pd.DataFrame,
    mcnemar_results: dict,
    timing_comparison: dict,
    output_path: pathlib.Path,
) -> None:
    """Generate a markdown comparison report."""
    with open(output_path, "w") as f:
        f.write("# km vs kmerdet Comparison Report\n\n")

        # Overall metrics
        f.write("## Overall Metrics\n\n")
        overall = summary_df[summary_df["category"] == "overall"]
        if not overall.empty:
            f.write("| Metric | km | kmerdet | Difference |\n")
            f.write("|--------|-----|---------|------------|\n")
            for _, row in overall.iterrows():
                km_val = f"{row['km']:.3f}" if row['km'] is not None else "N/A"
                kd_val = f"{row['kmerdet']:.3f}" if row['kmerdet'] is not None else "N/A"
                diff = ""
                if row['diff'] is not None:
                    d = row['diff']
                    arrow = "+" if d > 0.005 else ("-" if d < -0.005 else "~")
                    diff = f"{arrow}{abs(d):.3f}"
                f.write(f"| {row['metric']} | {km_val} | {kd_val} | {diff} |\n")
        f.write("\n")

        # Per-type breakdown
        f.write("## Sensitivity by Variant Type\n\n")
        per_type = summary_df[summary_df["category"].str.startswith("type_")]
        if not per_type.empty:
            f.write("| Variant Type | km | kmerdet | Difference |\n")
            f.write("|-------------|-----|---------|------------|\n")
            for _, row in per_type.iterrows():
                vtype = row["category"].replace("type_", "")
                km_val = f"{row['km']:.3f}" if row['km'] is not None else "N/A"
                kd_val = f"{row['kmerdet']:.3f}" if row['kmerdet'] is not None else "N/A"
                diff = ""
                if row['diff'] is not None:
                    d = row['diff']
                    arrow = "+" if d > 0.005 else ("-" if d < -0.005 else "~")
                    diff = f"{arrow}{abs(d):.3f}"
                f.write(f"| {vtype} | {km_val} | {kd_val} | {diff} |\n")
        f.write("\n")

        # Per-VAF breakdown
        f.write("## Sensitivity by VAF Bin\n\n")
        per_vaf = summary_df[summary_df["category"].str.startswith("vaf_")]
        if not per_vaf.empty:
            f.write("| VAF Bin | km | kmerdet | Difference |\n")
            f.write("|---------|-----|---------|------------|\n")
            for _, row in per_vaf.iterrows():
                label = row["category"].replace("vaf_", "")
                km_val = f"{row['km']:.3f}" if row['km'] is not None else "N/A"
                kd_val = f"{row['kmerdet']:.3f}" if row['kmerdet'] is not None else "N/A"
                diff = ""
                if row['diff'] is not None:
                    d = row['diff']
                    arrow = "+" if d > 0.005 else ("-" if d < -0.005 else "~")
                    diff = f"{arrow}{abs(d):.3f}"
                f.write(f"| {label} | {km_val} | {kd_val} | {diff} |\n")
        f.write("\n")

        # Concordance
        if not concordance_df.empty:
            f.write("## Concordance Matrix\n\n")
            counts = concordance_df["concordance"].value_counts()
            f.write("| Category | Count |\n")
            f.write("|----------|-------|\n")
            for cat in ["both_TP", "both_FN", "km_only", "kmerdet_only", "both_TN", "other"]:
                count = counts.get(cat, 0)
                f.write(f"| {cat} | {count} |\n")
            f.write("\n")
            total_present = counts.get("both_TP", 0) + counts.get("both_FN", 0) + \
                           counts.get("km_only", 0) + counts.get("kmerdet_only", 0)
            if total_present > 0:
                agreement = (counts.get("both_TP", 0) + counts.get("both_FN", 0)) / total_present
                f.write(f"**Agreement rate on present variants**: {agreement:.1%}\n\n")

        # McNemar's test
        if mcnemar_results:
            f.write("## Statistical Significance (McNemar's Test)\n\n")
            for suffix, result in mcnemar_results.items():
                f.write(f"### Parameters: {suffix}\n\n")
                f.write(f"- Variants detected by kmerdet only: {result.get('b_kmerdet_only', 0)}\n")
                f.write(f"- Variants detected by km only: {result.get('c_km_only', 0)}\n")
                f.write(f"- Both detected: {result.get('both_detected', 'N/A')}\n")
                f.write(f"- Both missed: {result.get('both_missed', 'N/A')}\n")
                f.write(f"- Chi-square statistic: {result.get('statistic', 0):.3f}\n")
                f.write(f"- p-value: {result.get('pvalue', 1.0):.4f}\n")
                f.write(f"- **{result.get('interpretation', '')}**\n\n")

        # Timing comparison
        if timing_comparison:
            f.write("## Performance\n\n")
            f.write("| Tool | Parameters | Wall Time (s) | Peak RSS (KB) |\n")
            f.write("|------|-----------|---------------|---------------|\n")
            for entry in timing_comparison:
                f.write(
                    f"| {entry['tool']} | {entry['params']} | "
                    f"{entry['wall_sec']:.2f} | "
                    f"{entry.get('peak_rss_kb', 'N/A')} |\n"
                )
        f.write("\n")

        f.write("---\n")
        f.write("*Generated by compare_results.py*\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--km-dir",      type=pathlib.Path, required=True,
                   help="Directory containing km benchmark outputs.")
    p.add_argument("--kmerdet-dir", type=pathlib.Path, required=True,
                   help="Directory containing kmerdet benchmark outputs.")
    p.add_argument("--truth",       type=pathlib.Path, required=True,
                   help="Ground truth TSV file.")
    p.add_argument("--output-dir",  type=pathlib.Path, default=pathlib.Path("comparison"),
                   help="Output directory [default: comparison/].")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover benchmark files
    km_files = discover_benchmark_files(args.km_dir)
    kmerdet_files = discover_benchmark_files(args.kmerdet_dir)

    if not km_files and not kmerdet_files:
        print("ERROR: No benchmark JSON files found in either directory.",
              file=sys.stderr)
        return 1

    print(f"Found {len(km_files)} km benchmark files, "
          f"{len(kmerdet_files)} kmerdet benchmark files.")

    # Build a map of suffix → (km_bench, kmerdet_bench)
    km_map = {suffix: path for suffix, path in km_files}
    kmerdet_map = {suffix: path for suffix, path in kmerdet_files}
    all_suffixes = sorted(set(km_map.keys()) | set(kmerdet_map.keys()))

    all_summaries = []
    all_concordance = pd.DataFrame()
    all_mcnemar = {}
    all_timing = []

    for suffix in all_suffixes:
        print(f"\n--- Comparing: {suffix} ---")

        km_bench = None
        kmerdet_bench = None
        km_metrics = {}
        kmerdet_metrics = {}
        km_per_type = pd.DataFrame()
        kmerdet_per_type = pd.DataFrame()
        km_per_vaf = pd.DataFrame()
        kmerdet_per_vaf = pd.DataFrame()

        if suffix in km_map:
            km_bench = load_benchmark_json(km_map[suffix])
            km_metrics = extract_metrics(km_bench)
            km_per_type = extract_per_type(km_bench)
            km_per_vaf = extract_per_vaf(km_bench)
            print(f"  km: sensitivity={km_metrics.get('sensitivity', 'N/A')}, "
                  f"TP={km_metrics.get('tp', 0)}, FP={km_metrics.get('fp', 0)}, "
                  f"FN={km_metrics.get('fn', 0)}")
        else:
            print("  km: no data")

        if suffix in kmerdet_map:
            kmerdet_bench = load_benchmark_json(kmerdet_map[suffix])
            kmerdet_metrics = extract_metrics(kmerdet_bench)
            kmerdet_per_type = extract_per_type(kmerdet_bench)
            kmerdet_per_vaf = extract_per_vaf(kmerdet_bench)
            print(f"  kmerdet: sensitivity={kmerdet_metrics.get('sensitivity', 'N/A')}, "
                  f"TP={kmerdet_metrics.get('tp', 0)}, FP={kmerdet_metrics.get('fp', 0)}, "
                  f"FN={kmerdet_metrics.get('fn', 0)}")
        else:
            print("  kmerdet: no data")

        # Build summary
        summary = generate_comparison_summary(
            km_metrics, kmerdet_metrics,
            km_per_type, kmerdet_per_type,
            km_per_vaf, kmerdet_per_vaf,
            suffix,
        )
        all_summaries.append(summary)

        # Build concordance from per-variant data
        if km_bench and kmerdet_bench:
            km_pv = extract_per_variant(km_bench)
            kmerdet_pv = extract_per_variant(kmerdet_bench)
            concordance = build_concordance(km_pv, kmerdet_pv)
            if not concordance.empty:
                concordance["params"] = suffix
                all_concordance = pd.concat([all_concordance, concordance],
                                            ignore_index=True)

            # McNemar's test on present variants
            if not km_pv.empty and not kmerdet_pv.empty:
                # Extract TP masks for present variants
                km_present = km_pv[km_pv.get("true_vaf", pd.Series(dtype=float)).astype(float) > 0]
                kd_present = kmerdet_pv[kmerdet_pv.get("true_vaf", pd.Series(dtype=float)).astype(float) > 0]

                if len(km_present) == len(kd_present) and len(km_present) > 0:
                    km_tp_mask = (km_present["status"] == "TP").to_numpy()
                    kd_tp_mask = (kd_present["status"] == "TP").to_numpy()
                    mcnemar = mcnemar_test(kd_tp_mask, km_tp_mask)
                    all_mcnemar[suffix] = mcnemar
                    print(f"  McNemar: b={mcnemar['b_kmerdet_only']}, "
                          f"c={mcnemar['c_km_only']}, "
                          f"p={mcnemar['pvalue']:.4f}")
                    print(f"  {mcnemar['interpretation']}")

        # Collect timing data
        for tool, tool_dir in [("km", args.km_dir), ("kmerdet", args.kmerdet_dir)]:
            timing_file = tool_dir / f"{tool}_timing_{suffix}.json"
            if timing_file.exists():
                with open(timing_file) as f:
                    t = json.load(f)
                all_timing.append({
                    "tool": tool,
                    "params": suffix,
                    "wall_sec": t.get("wall_sec", 0),
                    "peak_rss_kb": t.get("peak_rss_kb"),
                })

    # Write outputs
    if all_summaries:
        summary_df = pd.concat(all_summaries, ignore_index=True)
        summary_df.to_csv(out_dir / "comparison_summary.tsv", sep="\t", index=False)
        print(f"\nComparison summary: {out_dir / 'comparison_summary.tsv'}")
    else:
        summary_df = pd.DataFrame()

    if not all_concordance.empty:
        all_concordance.to_csv(out_dir / "concordance_matrix.tsv", sep="\t", index=False)
        print(f"Concordance matrix: {out_dir / 'concordance_matrix.tsv'}")

    # Save McNemar results
    if all_mcnemar:
        with open(out_dir / "mcnemar_results.json", "w") as f:
            json.dump(all_mcnemar, f, indent=2)

    # Generate report
    generate_report(
        summary_df=summary_df,
        concordance_df=all_concordance,
        mcnemar_results=all_mcnemar,
        timing_comparison=all_timing,
        output_path=out_dir / "comparison_report.md",
    )
    print(f"Comparison report: {out_dir / 'comparison_report.md'}")

    # Print final summary table
    print("\n=== km vs kmerdet Summary ===\n")
    if not summary_df.empty:
        overall = summary_df[summary_df["category"] == "overall"]
        if not overall.empty:
            print(f"{'Metric':<15} {'km':>8} {'kmerdet':>8} {'Diff':>8}")
            print("-" * 42)
            for _, row in overall.iterrows():
                km_v = f"{row['km']:.3f}" if row['km'] is not None else "N/A"
                kd_v = f"{row['kmerdet']:.3f}" if row['kmerdet'] is not None else "N/A"
                diff = f"{row['diff']:+.3f}" if row['diff'] is not None else ""
                print(f"{row['metric']:<15} {km_v:>8} {kd_v:>8} {diff:>8}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
