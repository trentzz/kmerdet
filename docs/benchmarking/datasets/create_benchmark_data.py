#!/usr/bin/env python3
"""create_benchmark_data.py — Generate jellyfish .jf benchmark databases.

Unlike generate_simulated.py (which produces mock TSVs), this script creates
actual jellyfish .jf database files by:
  1. Generating reference sequences with injected variants (SNV, ins, del, ITD)
  2. Tiling synthetic 150bp reads at controlled VAF and coverage
  3. Writing reads to FASTA, then running `jellyfish count` to produce .jf
  4. Writing target FASTA files and ground truth TSV

This allows head-to-head benchmarking of km and kmerdet on identical .jf inputs.

Presets:
  snv_indel                — 100+ mixed variants at 0.1%-50% VAF
  vaf_titration            — Fixed VAF levels (0.05% to 10%) for LOD characterization
  coverage_scaling         — Same variants at 500x, 1000x, 3000x, 5000x coverage
  large_indels             — 5-50bp INDELs at k=21 to stress INDEL sensitivity
  vaf_titration_standard   — Baseline sensitivity curve (0.05%-10% VAF, 160 variants)
  vaf_titration_ultra_low  — Ultra-low VAF frontier (0.0001%-0.1%, 280 variants)
  vaf_titration_extreme    — Extreme 0.0001% regime (100K/500Kx, 300 variants)
  coverage_vaf_matrix      — Coverage x VAF interaction (1K-100Kx, 600 variants)
  indel_size_sensitivity   — INDEL size vs detection (sizes 1-50bp, k=21/31)
  parameter_sensitivity    — Threshold optimization at fixed 0.05% VAF

Requirements:
  - jellyfish (>= 2.0) on $PATH
  - Python >= 3.9

Usage:
    python3 create_benchmark_data.py --preset snv_indel --output-dir /tmp/bench --seed 42
    python3 create_benchmark_data.py --preset all --base-dir /tmp/bench_all
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile
from typing import Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASES = "ACGT"
BASE_COMPLEMENT = str.maketrans("ACGT", "TGCA")
DEFAULT_READ_LENGTH = 150
DEFAULT_KMER_LENGTH = 31


# ---------------------------------------------------------------------------
# DNA utilities (shared with generate_simulated.py)
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    return seq.translate(BASE_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def random_dna(length: int, rng: random.Random, gc_fraction: float = 0.5) -> str:
    result = []
    for _ in range(length):
        r = rng.random()
        if r < gc_fraction / 2:
            result.append("G")
        elif r < gc_fraction:
            result.append("C")
        elif r < (1 + gc_fraction) / 2:
            result.append("A")
        else:
            result.append("T")
    return "".join(result)


def is_repetitive(seq: str, k: int, max_repeat_fraction: float = 0.3) -> bool:
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    if not kmers:
        return False
    unique = len(set(kmers))
    return (unique / len(kmers)) < (1 - max_repeat_fraction)


def homopolymer_run_length(seq: str) -> int:
    if not seq:
        return 0
    max_run = 1
    cur_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur_run += 1
            max_run = max(max_run, cur_run)
        else:
            cur_run = 1
    return max_run


def make_unique_target_sequence(
    length: int,
    kmer_length: int,
    rng: random.Random,
    max_homopolymer: int = 6,
    max_attempts: int = 100,
) -> str:
    for _ in range(max_attempts):
        seq = random_dna(length, rng, gc_fraction=0.45)
        if (not is_repetitive(seq, kmer_length)
                and homopolymer_run_length(seq) <= max_homopolymer):
            return seq
    return seq


# ---------------------------------------------------------------------------
# Variant injection (reused from generate_simulated.py)
# ---------------------------------------------------------------------------

def inject_snv(seq: str, pos: int, rng: random.Random) -> tuple[str, str, str]:
    ref_base = seq[pos]
    other_bases = [b for b in BASES if b != ref_base]
    alt_base = rng.choice(other_bases)
    alt_seq = seq[:pos] + alt_base + seq[pos + 1:]
    return alt_seq, ref_base, alt_base


def inject_insertion(seq: str, pos: int, insert_len: int, rng: random.Random) -> tuple[str, str, str]:
    anchor = seq[pos]
    inserted = random_dna(insert_len, rng)
    alt_seq = seq[:pos + 1] + inserted + seq[pos + 1:]
    return alt_seq, anchor, anchor + inserted


def inject_deletion(seq: str, pos: int, del_len: int) -> tuple[str, str, str]:
    anchor = seq[pos]
    deleted = seq[pos:pos + del_len + 1]
    alt_seq = seq[:pos + 1] + seq[pos + del_len + 1:]
    return alt_seq, deleted, anchor


def inject_itd(seq: str, pos: int, dup_len: int) -> tuple[str, str, str]:
    anchor = seq[pos]
    dup_seq = seq[pos + 1:pos + 1 + dup_len]
    if len(dup_seq) < dup_len:
        dup_seq += random_dna(dup_len - len(dup_seq), random.Random(pos))
    alt_seq = seq[:pos + 1] + dup_seq + seq[pos + 1:]
    return alt_seq, anchor, anchor + dup_seq


# ---------------------------------------------------------------------------
# Read tiling: generate synthetic FASTA reads
# ---------------------------------------------------------------------------

def tile_reads(
    seq: str,
    read_length: int,
    n_reads: int,
    rng: random.Random,
) -> list[str]:
    """Generate n_reads random-start reads from seq (both strands)."""
    reads = []
    seq_len = len(seq)
    if seq_len < read_length:
        # Sequence shorter than read length — return the whole thing
        reads.append(seq)
        return reads

    for _ in range(n_reads):
        start = rng.randint(0, seq_len - read_length)
        read = seq[start:start + read_length]
        # Randomly reverse-complement ~50% of reads
        if rng.random() < 0.5:
            read = reverse_complement(read)
        reads.append(read)
    return reads


def compute_n_reads(
    seq_len: int,
    read_length: int,
    coverage: float,
) -> int:
    """Number of reads to achieve target coverage over a sequence."""
    return max(1, int(math.ceil(coverage * seq_len / read_length)))


# ---------------------------------------------------------------------------
# Target generation with read tiling
# ---------------------------------------------------------------------------

def generate_target(
    target_id: str,
    variant_type: str,
    kmer_length: int,
    target_length: int,
    coverage: int,
    vaf: float,
    rng: random.Random,
    read_length: int = DEFAULT_READ_LENGTH,
    indel_len: int = 5,
    dup_len: int = 10,
) -> dict:
    """Generate a single target with variant and tiled reads."""
    ref_seq = make_unique_target_sequence(target_length, kmer_length, rng)

    # Choose variant position away from edges
    margin = kmer_length + 5
    if margin >= len(ref_seq) // 2:
        margin = max(5, len(ref_seq) // 4)
    pos = rng.randint(margin, len(ref_seq) - margin - max(indel_len, dup_len) - 2)

    # Inject variant
    if variant_type == "Substitution":
        alt_seq, ref_allele, alt_allele = inject_snv(ref_seq, pos, rng)
    elif variant_type == "Insertion":
        alt_seq, ref_allele, alt_allele = inject_insertion(ref_seq, pos, indel_len, rng)
    elif variant_type == "Deletion":
        del_end = min(pos + indel_len + 1, len(ref_seq))
        actual_del_len = del_end - pos - 1
        alt_seq, ref_allele, alt_allele = inject_deletion(ref_seq, pos, actual_del_len)
    elif variant_type == "ITD":
        dup_end = min(pos + dup_len + 1, len(ref_seq))
        actual_dup = dup_end - pos - 1
        alt_seq, ref_allele, alt_allele = inject_itd(ref_seq, pos, actual_dup)
    else:
        raise ValueError(f"Unknown variant type: {variant_type}")

    # Generate reads for ref and alt alleles
    ref_cov = coverage * (1.0 - vaf)
    alt_cov = coverage * vaf

    n_ref_reads = compute_n_reads(len(ref_seq), read_length, ref_cov)
    n_alt_reads = compute_n_reads(len(alt_seq), read_length, alt_cov) if vaf > 0 else 0

    ref_reads = tile_reads(ref_seq, read_length, n_ref_reads, rng)
    alt_reads = tile_reads(alt_seq, read_length, n_alt_reads, rng) if n_alt_reads > 0 else []

    return {
        "target_id":    target_id,
        "variant_type": variant_type,
        "ref_seq":      ref_seq,
        "alt_seq":      alt_seq,
        "ref_allele":   ref_allele,
        "alt_allele":   alt_allele,
        "pos":          pos + 1,       # 1-based
        "chrom":        f"sim{target_id}",
        "vaf":          vaf,
        "coverage":     coverage,
        "reads":        ref_reads + alt_reads,
        "n_ref_reads":  len(ref_reads),
        "n_alt_reads":  len(alt_reads),
    }


# ---------------------------------------------------------------------------
# Jellyfish database creation
# ---------------------------------------------------------------------------

def check_jellyfish() -> str:
    """Verify jellyfish is available and return its path."""
    jf_path = shutil.which("jellyfish")
    if jf_path is None:
        print("ERROR: 'jellyfish' not found on $PATH.", file=sys.stderr)
        print("Install jellyfish >= 2.0: https://github.com/gmarcais/Jellyfish",
              file=sys.stderr)
        sys.exit(1)

    # Verify version
    try:
        result = subprocess.run(
            [jf_path, "--version"],
            capture_output=True, text=True, timeout=10,
        )
        version_str = result.stdout.strip() or result.stderr.strip()
        print(f"  Using jellyfish: {jf_path} ({version_str})")
    except Exception as e:
        print(f"  Warning: could not get jellyfish version: {e}", file=sys.stderr)

    return jf_path


def create_jf_database(
    reads_fasta: pathlib.Path,
    output_jf: pathlib.Path,
    kmer_length: int,
    jellyfish_path: str,
    min_count: int = 2,
    hash_size: str = "100M",
) -> None:
    """Run jellyfish count to create a .jf database from a FASTA file."""
    cmd = [
        jellyfish_path, "count",
        "-m", str(kmer_length),
        "-C",                      # canonical mode (always)
        "-s", hash_size,
        "-L", str(min_count),      # minimum count threshold
        "-t", str(min(4, os.cpu_count() or 1)),
        "-o", str(output_jf),
        str(reads_fasta),
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"jellyfish count failed (exit {result.returncode})")
    print(f"  Created: {output_jf} ({output_jf.stat().st_size} bytes)")


# ---------------------------------------------------------------------------
# Dataset generation
# ---------------------------------------------------------------------------

def generate_dataset(
    output_dir: pathlib.Path,
    n_snvs: int,
    n_insertions: int,
    n_deletions: int,
    n_itds: int,
    n_absent: int,
    vaf_min: float,
    vaf_max: float,
    vaf_values: Optional[list[float]],
    coverage: int,
    kmer_length: int,
    target_length: int,
    read_length: int,
    indel_min_len: int,
    indel_max_len: int,
    dup_len: int,
    seed: int,
    jellyfish_path: str,
) -> None:
    """Generate a complete benchmark dataset with a real .jf database."""
    output_dir.mkdir(parents=True, exist_ok=True)
    targets_dir = output_dir / "targets"
    targets_dir.mkdir(exist_ok=True)

    rng = random.Random(seed)

    # Build variant list
    variant_specs: list[tuple[str, float]] = []
    total_real = n_snvs + n_insertions + n_deletions + n_itds

    if vaf_values:
        vaf_cycle = vaf_values * (total_real // len(vaf_values) + 1)
    else:
        vaf_cycle = None

    def next_vaf(i: int) -> float:
        if vaf_cycle:
            return vaf_cycle[i % len(vaf_cycle)]
        return math.exp(rng.uniform(math.log(vaf_min), math.log(vaf_max)))

    idx = 0
    for _ in range(n_snvs):
        variant_specs.append(("Substitution", next_vaf(idx)))
        idx += 1
    for _ in range(n_insertions):
        variant_specs.append(("Insertion", next_vaf(idx)))
        idx += 1
    for _ in range(n_deletions):
        variant_specs.append(("Deletion", next_vaf(idx)))
        idx += 1
    for _ in range(n_itds):
        variant_specs.append(("ITD", next_vaf(idx)))
        idx += 1
    for _ in range(n_absent):
        variant_specs.append(("Substitution", 0.0))
        idx += 1

    rng.shuffle(variant_specs)

    # Generate each target and collect all reads
    all_reads: list[str] = []
    ground_truth_rows: list[dict] = []
    fasta_records: list[tuple[str, str]] = []

    print(f"Generating {len(variant_specs)} targets in {output_dir}...")

    for i, (vtype, vaf) in enumerate(variant_specs):
        target_id = f"{i + 1:04d}"

        if vtype in ("Insertion", "Deletion"):
            indel_len = rng.randint(indel_min_len, indel_max_len)
        else:
            indel_len = 5

        dup = min(dup_len, target_length // 4)

        try:
            t = generate_target(
                target_id=target_id,
                variant_type=vtype if vaf > 0 else "Substitution",
                kmer_length=kmer_length,
                target_length=target_length,
                coverage=coverage,
                vaf=vaf,
                rng=rng,
                read_length=read_length,
                indel_len=indel_len,
                dup_len=dup,
            )
        except Exception as e:
            print(f"  Warning: could not generate target {target_id}: {e}",
                  file=sys.stderr)
            continue

        # Write target FASTA (reference sequence only — what km/kmerdet detects against)
        fasta_path = targets_dir / f"target_{target_id}.fa"
        with open(fasta_path, "w") as f:
            f.write(f">{t['target_id']}\n{t['ref_seq']}\n")
        fasta_records.append((t["target_id"], t["ref_seq"]))

        # Collect reads for the jellyfish database
        all_reads.extend(t["reads"])

        # Ground truth
        ground_truth_rows.append({
            "chrom":        t["chrom"],
            "pos":          t["pos"],
            "ref_allele":   t["ref_allele"],
            "alt_allele":   t["alt_allele"],
            "variant_type": vtype,
            "true_vaf":     vaf,
            "category":     vtype,
        })

    # Write combined targets FASTA
    all_targets_path = output_dir / "all_targets.fa"
    with open(all_targets_path, "w") as f:
        for tid, seq in fasta_records:
            f.write(f">{tid}\n{seq}\n")

    # Write ground truth TSV
    truth_path = output_dir / "ground_truth.tsv"
    print(f"Writing {len(ground_truth_rows)} truth entries to ground_truth.tsv...")
    with open(truth_path, "w") as f:
        f.write("chrom\tpos\tref\talt\ttype\ttrue_vaf\tcategory\n")
        for row in ground_truth_rows:
            f.write(
                f"{row['chrom']}\t{row['pos']}\t{row['ref_allele']}\t"
                f"{row['alt_allele']}\t{row['variant_type']}\t"
                f"{row['true_vaf']:.6f}\t{row['category']}\n"
            )

    # Write reads to temp FASTA and run jellyfish
    print(f"Writing {len(all_reads)} reads to temporary FASTA...")
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False, dir=str(output_dir)
    ) as tmp_reads:
        for i, read in enumerate(all_reads):
            tmp_reads.write(f">read_{i}\n{read}\n")
        tmp_reads_path = pathlib.Path(tmp_reads.name)

    try:
        jf_path = output_dir / "database.jf"
        create_jf_database(
            reads_fasta=tmp_reads_path,
            output_jf=jf_path,
            kmer_length=kmer_length,
            jellyfish_path=jellyfish_path,
            min_count=2,
        )
    finally:
        # Clean up temp reads file
        tmp_reads_path.unlink(missing_ok=True)

    # Write manifest
    n_present = sum(1 for r in ground_truth_rows if r["true_vaf"] > 0)
    manifest = {
        "generated_by": "create_benchmark_data.py",
        "seed": seed,
        "kmer_length": kmer_length,
        "target_length": target_length,
        "read_length": read_length,
        "coverage": coverage,
        "vaf_min": vaf_min,
        "vaf_max": vaf_max,
        "vaf_values": vaf_values,
        "n_snvs": n_snvs,
        "n_insertions": n_insertions,
        "n_deletions": n_deletions,
        "n_itds": n_itds,
        "n_absent": n_absent,
        "indel_min_len": indel_min_len,
        "indel_max_len": indel_max_len,
        "dup_len": dup_len,
        "total_targets": len(fasta_records),
        "total_truth_entries": len(ground_truth_rows),
        "n_present": n_present,
        "n_absent_truth": len(ground_truth_rows) - n_present,
        "total_reads": len(all_reads),
        "database_type": "jellyfish",
        "files": {
            "targets_dir": str(targets_dir),
            "all_targets_fasta": str(all_targets_path),
            "database_jf": str(jf_path),
            "ground_truth": str(truth_path),
        },
    }
    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nDataset generated in: {output_dir}")
    print(f"  Targets:     {len(fasta_records)}")
    print(f"  Truth variants (present): {n_present}")
    print(f"  Truth variants (absent):  {len(ground_truth_rows) - n_present}")
    print(f"  Total reads: {len(all_reads)}")
    print(f"  Database:    {jf_path}")
    print(f"  Files:")
    for name, path in manifest["files"].items():
        print(f"    {name}: {path}")
    print()


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

PRESETS = {
    "snv_indel": dict(
        description="Standard SNV + INDEL accuracy dataset (100+ variants)",
        n_snvs=50,
        n_insertions=25,
        n_deletions=25,
        n_itds=10,
        n_absent=15,
        vaf_min=0.001,
        vaf_max=0.5,
        coverage=3000,
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=10,
    ),
    "vaf_titration": dict(
        description="Fixed VAF levels for LOD characterization",
        n_snvs=20,
        n_insertions=10,
        n_deletions=10,
        n_itds=5,
        n_absent=10,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10],
        coverage=3000,
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "coverage_scaling": dict(
        description="Same variants at multiple coverage depths",
        n_snvs=20,
        n_insertions=10,
        n_deletions=10,
        n_itds=5,
        n_absent=5,
        vaf_min=0.005,
        vaf_max=0.1,
        coverage=3000,  # overridden per-run below
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=10,
    ),
    "large_indels": dict(
        description="Large INDELs (5-50bp) at k=21",
        n_snvs=0,
        n_insertions=20,
        n_deletions=20,
        n_itds=10,
        n_absent=5,
        vaf_min=0.01,
        vaf_max=0.2,
        coverage=3000,
        kmer_length=21,
        target_length=200,
        read_length=150,
        indel_min_len=5,
        indel_max_len=50,
    ),
    "vaf_titration_standard": dict(
        description="Baseline sensitivity curve (0.05%-10% VAF, 3000x, 160 variants)",
        n_snvs=10,
        n_insertions=5,
        n_deletions=5,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005],
        coverage=3000,
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "vaf_titration_ultra_low": dict(
        description="Ultra-low VAF frontier (0.0001%-0.1% VAF, 50000x, 280 variants)",
        n_snvs=20,
        n_insertions=10,
        n_deletions=10,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001],
        coverage=50000,
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "vaf_titration_extreme": dict(
        description="Extreme VAF regime (0.0001%, 100000x/500000x, 300 variants)",
        n_snvs=50,
        n_insertions=25,
        n_deletions=25,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.00001, 0.000005, 0.000001],
        coverage=100000,  # overridden per-run via EXTREME_DEPTHS
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "coverage_vaf_matrix": dict(
        description="Coverage x VAF interaction matrix (1K-100Kx, 600 variants)",
        n_snvs=20,
        n_insertions=10,
        n_deletions=10,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.001, 0.0001, 0.00001, 0.000001],
        coverage=3000,  # overridden per-run via COVERAGE_VAF_MATRIX_DEPTHS
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "indel_size_sensitivity": dict(
        description="INDEL size vs detection sensitivity (sizes 1-50bp, k=21/31, 800 variants)",
        n_snvs=0,
        n_insertions=10,
        n_deletions=10,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.01, 0.001],
        coverage=3000,
        kmer_length=21,
        target_length=200,
        read_length=150,
        indel_min_len=1,
        indel_max_len=50,
    ),
    "parameter_sensitivity": dict(
        description="Threshold optimization at fixed 0.05% VAF (10000x, 100 variants)",
        n_snvs=50,
        n_insertions=25,
        n_deletions=25,
        n_itds=0,
        n_absent=0,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.0005],
        coverage=10000,
        kmer_length=31,
        target_length=180,
        read_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
}

# Coverage scaling preset generates multiple datasets at different depths
COVERAGE_SCALING_DEPTHS = [500, 1000, 3000, 5000]

# Coverage × VAF matrix generates one .jf per coverage level
COVERAGE_VAF_MATRIX_DEPTHS = [1000, 3000, 10000, 50000, 100000]

# Extreme VAF titration generates at these ultra-high depths
EXTREME_DEPTHS = [100000, 500000]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "--output-dir", "-o",
        type=pathlib.Path,
        help="Output directory for the generated dataset. "
             "Required unless --preset all is used with --base-dir.",
    )
    p.add_argument(
        "--base-dir",
        type=pathlib.Path,
        default=pathlib.Path("data/benchmark"),
        help="Base directory for --preset all. Each preset gets a subdirectory. "
             "[default: data/benchmark/]",
    )
    p.add_argument(
        "--preset",
        choices=list(PRESETS.keys()) + ["all"],
        help="Use a predefined dataset configuration.",
    )

    # Variant counts
    p.add_argument("--n-snvs",       type=int, default=50,  help="Number of SNVs [default: 50]")
    p.add_argument("--n-insertions", type=int, default=25,  help="Number of insertions [default: 25]")
    p.add_argument("--n-deletions",  type=int, default=25,  help="Number of deletions [default: 25]")
    p.add_argument("--n-itds",       type=int, default=10,  help="Number of ITDs [default: 10]")
    p.add_argument("--n-absent",     type=int, default=15,
                   help="Number of absent variants (specificity controls) [default: 15]")

    # VAF settings
    p.add_argument("--vaf-min",    type=float, default=0.001,
                   help="Minimum VAF for log-uniform sampling [default: 0.001]")
    p.add_argument("--vaf-max",    type=float, default=0.5,
                   help="Maximum VAF for log-uniform sampling [default: 0.5]")
    p.add_argument("--vaf-values", type=str,
                   help="Explicit comma-separated VAF values (overrides --vaf-min/max).")

    # Sequencing parameters
    p.add_argument("--coverage",      type=int, default=3000,
                   help="Simulated sequencing depth [default: 3000]")
    p.add_argument("--kmer-length",   type=int, default=31,
                   help="K-mer length [default: 31]")
    p.add_argument("--target-length", type=int, default=180,
                   help="Target sequence length (bp) [default: 180]")
    p.add_argument("--read-length",   type=int, default=150,
                   help="Synthetic read length (bp) [default: 150]")

    # INDEL settings
    p.add_argument("--indel-min-len", type=int, default=1,
                   help="Minimum INDEL length (bp) [default: 1]")
    p.add_argument("--indel-max-len", type=int, default=10,
                   help="Maximum INDEL length (bp) [default: 10]")
    p.add_argument("--dup-len",       type=int, default=10,
                   help="ITD duplication length (bp) [default: 10]")

    p.add_argument("--seed", type=int, default=42, help="Random seed [default: 42]")

    return p


def run_preset(
    name: str,
    cfg: dict,
    output_dir: pathlib.Path,
    seed: int,
    jellyfish_path: str,
    vaf_values_override: Optional[list[float]] = None,
) -> None:
    """Run a single preset configuration."""
    vaf_vals = cfg.get("vaf_values") or vaf_values_override
    vaf_min = cfg.get("vaf_min") or 0.001
    vaf_max = cfg.get("vaf_max") or 0.5

    generate_dataset(
        output_dir=output_dir,
        n_snvs=cfg["n_snvs"],
        n_insertions=cfg["n_insertions"],
        n_deletions=cfg["n_deletions"],
        n_itds=cfg["n_itds"],
        n_absent=cfg["n_absent"],
        vaf_min=vaf_min,
        vaf_max=vaf_max,
        vaf_values=vaf_vals,
        coverage=cfg["coverage"],
        kmer_length=cfg["kmer_length"],
        target_length=cfg["target_length"],
        read_length=cfg["read_length"],
        indel_min_len=cfg["indel_min_len"],
        indel_max_len=cfg["indel_max_len"],
        dup_len=cfg.get("dup_len", 10),
        seed=seed,
        jellyfish_path=jellyfish_path,
    )


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # Check jellyfish availability
    jellyfish_path = check_jellyfish()

    # Parse explicit VAF values if provided
    vaf_values: Optional[list[float]] = None
    if args.vaf_values:
        vaf_values = [float(v.strip()) for v in args.vaf_values.split(",")]

    if args.preset == "all":
        base = args.base_dir
        for name, cfg in PRESETS.items():
            print(f"\n{'='*60}")
            print(f"=== Generating preset: {name} ===")
            print(f"    {cfg['description']}")
            print(f"{'='*60}")

            if name == "coverage_scaling":
                # Generate multiple datasets at different coverage depths
                for depth in COVERAGE_SCALING_DEPTHS:
                    out = base / f"{name}_{depth}x"
                    cfg_copy = dict(cfg)
                    cfg_copy["coverage"] = depth
                    run_preset(name, cfg_copy, out, args.seed, jellyfish_path,
                               vaf_values)
            elif name == "coverage_vaf_matrix":
                # Generate one dataset per coverage level
                for depth in COVERAGE_VAF_MATRIX_DEPTHS:
                    out = base / f"{name}_{depth}x"
                    cfg_copy = dict(cfg)
                    cfg_copy["coverage"] = depth
                    run_preset(name, cfg_copy, out, args.seed, jellyfish_path,
                               vaf_values)
            elif name == "vaf_titration_extreme":
                # Generate at each extreme depth
                for depth in EXTREME_DEPTHS:
                    out = base / f"{name}_{depth}x"
                    cfg_copy = dict(cfg)
                    cfg_copy["coverage"] = depth
                    run_preset(name, cfg_copy, out, args.seed, jellyfish_path,
                               vaf_values)
            elif name == "indel_size_sensitivity":
                # Generate sub-datasets at different INDEL sizes and k values
                indel_sizes = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]
                k_values = [21, 31]
                for indel_size in indel_sizes:
                    for k in k_values:
                        out = base / f"{name}_size{indel_size}_k{k}"
                        cfg_copy = dict(cfg)
                        cfg_copy["indel_min_len"] = indel_size
                        cfg_copy["indel_max_len"] = indel_size
                        cfg_copy["kmer_length"] = k
                        run_preset(name, cfg_copy, out, args.seed,
                                   jellyfish_path, vaf_values)
            else:
                out = base / name
                run_preset(name, cfg, out, args.seed, jellyfish_path, vaf_values)
        return 0

    elif args.preset and args.preset in PRESETS:
        cfg = PRESETS[args.preset]
        out = args.output_dir or (args.base_dir / args.preset)

        if args.preset == "coverage_scaling":
            for depth in COVERAGE_SCALING_DEPTHS:
                out_depth = out.parent / f"{out.name}_{depth}x"
                cfg_copy = dict(cfg)
                cfg_copy["coverage"] = depth
                run_preset(args.preset, cfg_copy, out_depth, args.seed,
                           jellyfish_path, vaf_values)
        elif args.preset == "coverage_vaf_matrix":
            for depth in COVERAGE_VAF_MATRIX_DEPTHS:
                out_depth = out.parent / f"{out.name}_{depth}x"
                cfg_copy = dict(cfg)
                cfg_copy["coverage"] = depth
                run_preset(args.preset, cfg_copy, out_depth, args.seed,
                           jellyfish_path, vaf_values)
        elif args.preset == "vaf_titration_extreme":
            for depth in EXTREME_DEPTHS:
                out_depth = out.parent / f"{out.name}_{depth}x"
                cfg_copy = dict(cfg)
                cfg_copy["coverage"] = depth
                run_preset(args.preset, cfg_copy, out_depth, args.seed,
                           jellyfish_path, vaf_values)
        elif args.preset == "indel_size_sensitivity":
            indel_sizes = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]
            k_values = [21, 31]
            for indel_size in indel_sizes:
                for k in k_values:
                    out_sub = out.parent / f"{out.name}_size{indel_size}_k{k}"
                    cfg_copy = dict(cfg)
                    cfg_copy["indel_min_len"] = indel_size
                    cfg_copy["indel_max_len"] = indel_size
                    cfg_copy["kmer_length"] = k
                    run_preset(args.preset, cfg_copy, out_sub, args.seed,
                               jellyfish_path, vaf_values)
        else:
            run_preset(args.preset, cfg, out, args.seed, jellyfish_path,
                       vaf_values)
        return 0

    else:
        # Custom configuration
        if args.output_dir is None:
            print("ERROR: --output-dir is required (or use --preset)",
                  file=sys.stderr)
            return 1
        generate_dataset(
            output_dir=args.output_dir,
            n_snvs=args.n_snvs,
            n_insertions=args.n_insertions,
            n_deletions=args.n_deletions,
            n_itds=args.n_itds,
            n_absent=args.n_absent,
            vaf_min=args.vaf_min,
            vaf_max=args.vaf_max,
            vaf_values=vaf_values,
            coverage=args.coverage,
            kmer_length=args.kmer_length,
            target_length=args.target_length,
            read_length=args.read_length,
            indel_min_len=args.indel_min_len,
            indel_max_len=args.indel_max_len,
            dup_len=args.dup_len,
            seed=args.seed,
            jellyfish_path=jellyfish_path,
        )
        return 0


if __name__ == "__main__":
    sys.exit(main())
