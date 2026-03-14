#!/usr/bin/env python3
"""parse_km_output.py — Normalize km find_mutation output to kmerdet-compatible format.

km find_mutation writes a 12-column TSV to stdout (defined by PathQuant.__str__()
in km's source). This script converts it to the 19-column kmerdet TSV format so
that both tools' outputs can be compared with a single analysis pipeline.

km output columns (from PathQuant.__str__()):
  1. Database       — jellyfish database path (→ sample)
  2. Query          — target FASTA path (→ target)
  3. Type           — variant type (Reference, Substitution, Insertion, Deletion, ITD)
  4. Variant_name   — e.g. "41:A/T:41" for SNV at pos 41 (ref A → alt T)
  5. rVAF           — relative VAF (NaN for Reference calls)
  6. Expression     — expression estimate (NNLS coefficient)
  7. Min_coverage   — minimum k-mer count on path
  8. Start_kmer_count — count of the starting k-mer
  9. Ref_sequence   — reference path DNA sequence
  10. Alt_sequence  — alternate path DNA sequence
  11. Info          — additional annotations

kmerdet TSV columns:
  1-11: same semantics as km (sample, target, variant_type, variant_name,
        rvaf, expression, min_coverage, start_kmer_count,
        ref_sequence, alt_sequence, info)
  12-15: chrom, pos, ref_allele, alt_allele (parsed from variant_name)
  16-19: pvalue, qual, ci_lower, ci_upper (set to NA for km output)

Usage:
    python3 parse_km_output.py --input km_output.tsv --output km_normalized.tsv
    km find_mutation -t targets/ db.jf | python3 parse_km_output.py --output km_normalized.tsv
"""

from __future__ import annotations

import argparse
import math
import os
import pathlib
import re
import sys
from typing import Optional


# ---------------------------------------------------------------------------
# km variant_name parsing
# ---------------------------------------------------------------------------

def parse_variant_name(variant_name: str, variant_type: str) -> dict:
    """Parse km's variant_name field to extract position, ref, and alt alleles.

    km variant_name formats (from diff_path_without_overlap):
      SNV:       "41:A/T:41"        → pos=41, ref=A, alt=T
      Insertion: "41:ins20:61"       → pos=41, ref=., alt=(from alt_sequence)
      Deletion:  "41:del20:61"       → pos=41, ref=(from ref_sequence), alt=.
      ITD:       "41:ITD20:61"       → pos=41, ref=., alt=(duplicated seq)
      Complex:   "41:A/T|42:ins3:45" → first component
      Reference: ""                  → no variant

    Returns dict with keys: pos, ref_allele, alt_allele (may be None for complex).
    """
    result = {"pos": None, "ref_allele": None, "alt_allele": None}

    if not variant_name or variant_name.strip() == "" or variant_type == "Reference":
        return result

    # Handle complex variants (pipe-separated): take first component
    if "|" in variant_name:
        variant_name = variant_name.split("|")[0]

    # SNV format: "pos:ref/alt:pos"
    snv_match = re.match(r"^(\d+):([ACGT]+)/([ACGT]+):(\d+)$", variant_name)
    if snv_match:
        result["pos"] = int(snv_match.group(1))
        result["ref_allele"] = snv_match.group(2)
        result["alt_allele"] = snv_match.group(3)
        return result

    # Insertion format: "pos:insN:pos"
    ins_match = re.match(r"^(\d+):ins(\d+):(\d+)$", variant_name)
    if ins_match:
        result["pos"] = int(ins_match.group(1))
        result["ref_allele"] = "."
        result["alt_allele"] = f"ins{ins_match.group(2)}"
        return result

    # Deletion format: "pos:delN:pos"
    del_match = re.match(r"^(\d+):del(\d+):(\d+)$", variant_name)
    if del_match:
        result["pos"] = int(del_match.group(1))
        result["ref_allele"] = f"del{del_match.group(2)}"
        result["alt_allele"] = "."
        return result

    # ITD format: "pos:ITDN:pos"
    itd_match = re.match(r"^(\d+):ITD(\d+):(\d+)$", variant_name)
    if itd_match:
        result["pos"] = int(itd_match.group(1))
        result["ref_allele"] = "."
        result["alt_allele"] = f"ITD{itd_match.group(2)}"
        return result

    # Fallback: try to parse just the leading position
    pos_match = re.match(r"^(\d+):", variant_name)
    if pos_match:
        result["pos"] = int(pos_match.group(1))

    return result


def extract_sample_name(db_path: str) -> str:
    """Extract a sample name from the jellyfish database path.

    km writes the full path to the database file. We extract just the
    filename stem (e.g., "/path/to/patient01.jf" → "patient01").
    """
    return pathlib.Path(db_path).stem


def extract_target_name(query_path: str) -> str:
    """Extract target name from the target FASTA path.

    km writes the full path. We extract the filename stem.
    """
    return pathlib.Path(query_path).stem


def normalize_rvaf(rvaf_str: str) -> str:
    """Normalize rVAF: km writes 'nan' for Reference calls → '0.0'."""
    try:
        val = float(rvaf_str)
        if math.isnan(val):
            return "0.0"
        return f"{val:.6f}"
    except (ValueError, TypeError):
        return "0.0"


# ---------------------------------------------------------------------------
# Main conversion logic
# ---------------------------------------------------------------------------

def convert_km_line(line: str) -> Optional[str]:
    """Convert a single km output line to kmerdet TSV format.

    Returns None for comment lines or unparseable lines.
    """
    line = line.rstrip("\n\r")

    # Skip comment lines (km writes # comments sometimes)
    if line.startswith("#"):
        return None

    # Skip empty lines
    if not line.strip():
        return None

    fields = line.split("\t")
    if len(fields) < 11:
        print(f"  Warning: skipping short line ({len(fields)} fields): {line[:80]}...",
              file=sys.stderr)
        return None

    # km columns (0-indexed)
    database       = fields[0]
    query          = fields[1]
    variant_type   = fields[2]
    variant_name   = fields[3]
    rvaf           = fields[4]
    expression     = fields[5]
    min_coverage   = fields[6]
    start_kmer_cnt = fields[7]
    ref_sequence   = fields[8]
    alt_sequence   = fields[9]
    info           = fields[10] if len(fields) > 10 else ""

    # Normalize fields
    sample = extract_sample_name(database)
    target = extract_target_name(query)
    rvaf_norm = normalize_rvaf(rvaf)

    # Parse variant_name for position and alleles
    parsed = parse_variant_name(variant_name, variant_type)
    chrom = target  # For simulated data, target name serves as chrom
    pos = str(parsed["pos"]) if parsed["pos"] is not None else "NA"
    ref_allele = parsed["ref_allele"] or "NA"
    alt_allele = parsed["alt_allele"] or "NA"

    # For INDELs where we got placeholder alleles, try to extract from sequences
    if ref_allele.startswith("del") and ref_sequence and alt_sequence:
        # Extract deleted bases from sequence difference
        pass  # Keep the symbolic representation for matching
    if alt_allele.startswith("ins") and ref_sequence and alt_sequence:
        pass  # Keep the symbolic representation

    # kmerdet extra columns (not available in km output)
    pvalue = "NA"
    qual = "NA"
    ci_lower = "NA"
    ci_upper = "NA"

    # Build 19-column kmerdet line
    out_fields = [
        sample, target, variant_type, variant_name,
        rvaf_norm, expression, min_coverage, start_kmer_cnt,
        ref_sequence, alt_sequence, info,
        chrom, pos, ref_allele, alt_allele,
        pvalue, qual, ci_lower, ci_upper,
    ]
    return "\t".join(out_fields)


def convert_km_output(
    input_lines: list[str],
    skip_header: bool = True,
) -> list[str]:
    """Convert all km output lines to kmerdet TSV format."""
    output_lines = []

    # Write kmerdet header
    header = "\t".join([
        "sample", "target", "variant_type", "variant_name",
        "rvaf", "expression", "min_coverage", "start_kmer_count",
        "ref_sequence", "alt_sequence", "info",
        "chrom", "pos", "ref_allele", "alt_allele",
        "pvalue", "qual", "ci_lower", "ci_upper",
    ])
    output_lines.append(header)

    for i, line in enumerate(input_lines):
        # Skip km header (first non-comment line if it looks like a header)
        if skip_header and i == 0:
            lower = line.lower()
            if "database" in lower or "query" in lower or "type" in lower:
                continue

        converted = convert_km_line(line)
        if converted is not None:
            output_lines.append(converted)

    return output_lines


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input", "-i",
        type=pathlib.Path,
        help="km find_mutation output TSV file. Reads from stdin if omitted.",
    )
    p.add_argument(
        "--output", "-o",
        type=pathlib.Path,
        help="Output normalized TSV file. Writes to stdout if omitted.",
    )
    p.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write a header line to output.",
    )
    p.add_argument(
        "--sample-name",
        type=str,
        help="Override sample name (instead of extracting from Database column).",
    )
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # Read input
    if args.input:
        if not args.input.exists():
            print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
            return 1
        with open(args.input) as f:
            lines = f.readlines()
    else:
        lines = sys.stdin.readlines()

    if not lines:
        print("WARNING: Empty input.", file=sys.stderr)
        return 0

    # Convert
    output_lines = convert_km_output(lines)

    if args.no_header and output_lines:
        output_lines = output_lines[1:]  # Remove header

    # Apply sample name override
    if args.sample_name and output_lines:
        new_lines = [output_lines[0]]  # Keep header
        for line in output_lines[1:]:
            fields = line.split("\t")
            if fields:
                fields[0] = args.sample_name
            new_lines.append("\t".join(fields))
        output_lines = new_lines

    # Write output
    output_text = "\n".join(output_lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(output_text)
        print(f"Normalized output written to: {args.output}", file=sys.stderr)
        print(f"  Lines: {len(output_lines) - 1} (excluding header)", file=sys.stderr)
    else:
        sys.stdout.write(output_text)

    return 0


if __name__ == "__main__":
    sys.exit(main())
