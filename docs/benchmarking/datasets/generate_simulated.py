#!/usr/bin/env python3
"""generate_simulated.py — Generate synthetic benchmark datasets for kmerdet.

Creates:
  - targets/          FASTA files (one per target) with reference sequences
  - mock_db.tsv       k-mer count table (tab-separated: kmer<TAB>count)
  - ground_truth.tsv  Ground truth variants with true VAF
  - manifest.json     Run parameters for reproducibility

The mock database can be used directly with kmerdet when jellyfish is not
installed (kmerdet will detect the .tsv extension and use its mock reader).

Usage examples:
    # Standard SNV + INDEL dataset
    python3 generate_simulated.py --output-dir /tmp/bench/snv_indel

    # Ultra-low-VAF dataset (0.01% - 0.1%)
    python3 generate_simulated.py \\
        --output-dir /tmp/bench/ultra_low_vaf \\
        --vaf-min 0.0001 --vaf-max 0.001 \\
        --coverage 5000

    # Large INDEL dataset (stresses k-mer length constraint)
    python3 generate_simulated.py \\
        --output-dir /tmp/bench/large_indels \\
        --n-snvs 0 --n-insertions 10 --n-deletions 10 \\
        --indel-min-len 5 --indel-max-len 50 \\
        --kmer-length 21

    # Full preset suite
    python3 generate_simulated.py --preset all --base-dir /tmp/bench
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import random
import sys
from typing import Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASES = "ACGT"
BASE_COMPLEMENT = str.maketrans("ACGT", "TGCA")


# ---------------------------------------------------------------------------
# DNA utilities
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    return seq.translate(BASE_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    """Return the lexicographically smaller of kmer and its reverse complement."""
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def random_dna(length: int, rng: random.Random, gc_fraction: float = 0.5) -> str:
    """Generate a random DNA sequence with approximate GC content."""
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
    """Return True if sequence has too many repeated k-mers (ambiguous anchors)."""
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    if not kmers:
        return False
    unique = len(set(kmers))
    return (unique / len(kmers)) < (1 - max_repeat_fraction)


def homopolymer_run_length(seq: str) -> int:
    """Return the length of the longest homopolymer run in seq."""
    if not seq:
        return 0
    max_run = 1
    cur_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
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
    """Generate a random target sequence that is not too repetitive or homopolymeric."""
    for _ in range(max_attempts):
        seq = random_dna(length, rng, gc_fraction=0.45)
        if (not is_repetitive(seq, kmer_length) and
                homopolymer_run_length(seq) <= max_homopolymer):
            return seq
    # Give up and return last attempt (may be suboptimal)
    return seq


# ---------------------------------------------------------------------------
# Variant creation
# ---------------------------------------------------------------------------

def inject_snv(seq: str, pos: int, rng: random.Random) -> tuple[str, str, str]:
    """Replace base at pos with a different base. Returns (alt_seq, ref_allele, alt_allele)."""
    ref_base = seq[pos]
    other_bases = [b for b in BASES if b != ref_base]
    alt_base = rng.choice(other_bases)
    alt_seq = seq[:pos] + alt_base + seq[pos+1:]
    return alt_seq, ref_base, alt_base


def inject_insertion(seq: str, pos: int, insert_len: int, rng: random.Random) -> tuple[str, str, str]:
    """Insert random sequence after pos. Returns (alt_seq, ref_allele, alt_allele).

    VCF-style: ref = seq[pos], alt = seq[pos] + inserted_sequence
    """
    anchor = seq[pos]
    inserted = random_dna(insert_len, rng)
    alt_seq = seq[:pos+1] + inserted + seq[pos+1:]
    return alt_seq, anchor, anchor + inserted


def inject_deletion(seq: str, pos: int, del_len: int) -> tuple[str, str, str]:
    """Delete del_len bases starting at pos+1. Returns (alt_seq, ref_allele, alt_allele).

    VCF-style: ref = seq[pos:pos+del_len+1], alt = seq[pos]
    """
    anchor = seq[pos]
    deleted = seq[pos:pos+del_len+1]   # anchor + deleted bases
    alt_seq = seq[:pos+1] + seq[pos+del_len+1:]
    return alt_seq, deleted, anchor


def inject_itd(seq: str, pos: int, dup_len: int) -> tuple[str, str, str]:
    """Tandem duplication of dup_len bases starting at pos+1.

    The inserted sequence is a copy of seq[pos+1:pos+1+dup_len].
    VCF-style: ref = seq[pos], alt = seq[pos] + duplicated_sequence
    """
    anchor = seq[pos]
    dup_seq = seq[pos+1:pos+1+dup_len]
    if len(dup_seq) < dup_len:
        # Not enough sequence; use random
        dup_seq += random_dna(dup_len - len(dup_seq), random.Random(pos))
    alt_seq = seq[:pos+1] + dup_seq + seq[pos+1:]
    return alt_seq, anchor, anchor + dup_seq


# ---------------------------------------------------------------------------
# K-mer counting (simulated counts from Poisson model)
# ---------------------------------------------------------------------------

def simulate_kmer_counts(
    ref_seq: str,
    alt_seq: str,
    kmer_length: int,
    coverage: int,
    vaf: float,
    rng: random.Random,
    error_rate: float = 0.001,
) -> dict[str, int]:
    """Simulate k-mer counts for a target region with a known variant.

    Models:
    - Counts from reference allele reads: Poisson(coverage * (1 - vaf))
    - Counts from alt allele reads: Poisson(coverage * vaf)
    - Canonical counting: each k-mer is stored as min(kmer, revcomp(kmer))
    - Error k-mers: each base independently mutated at error_rate, added as singletons

    Returns a dict: canonical_kmer -> count.
    """
    counts: dict[str, int] = {}

    def add_kmer_counts(seq: str, mean_cov: float) -> None:
        n = len(seq)
        for i in range(n - kmer_length + 1):
            kmer = seq[i:i+kmer_length]
            if len(kmer) < kmer_length:
                continue
            ck = canonical(kmer)
            c = rng.randint(0, 2 * int(mean_cov)) if mean_cov > 0 else 0  # simplified Poisson
            if c > 0:
                counts[ck] = counts.get(ck, 0) + c

    # Reference k-mers
    add_kmer_counts(ref_seq, coverage * (1.0 - vaf))

    # Variant k-mers
    if vaf > 0:
        add_kmer_counts(alt_seq, coverage * vaf)

    # Sequencing error k-mers (very low count, mostly singletons)
    for seq in (ref_seq, alt_seq):
        for i in range(len(seq) - kmer_length + 1):
            kmer = seq[i:i+kmer_length]
            for pos in range(len(kmer)):
                if rng.random() < error_rate:
                    other = rng.choice([b for b in BASES if b != kmer[pos]])
                    err_kmer = kmer[:pos] + other + kmer[pos+1:]
                    ck = canonical(err_kmer)
                    if ck not in counts:
                        # Discard singleton error k-mers (jellyfish -L 2 would remove these)
                        pass  # intentionally do not add count-1 error k-mers by default

    # Remove k-mers with count < 2 (mimics jellyfish -L 2)
    return {k: v for k, v in counts.items() if v >= 2}


# ---------------------------------------------------------------------------
# Target generation
# ---------------------------------------------------------------------------

def generate_target(
    target_id: str,
    variant_type: str,
    kmer_length: int,
    target_length: int,
    coverage: int,
    vaf: float,
    rng: random.Random,
    indel_len: int = 5,
    dup_len: int = 10,
) -> dict:
    """Generate a single synthetic target with one injected variant.

    Returns a dict with:
      ref_seq, alt_seq, variant_type, ref_allele, alt_allele,
      pos (1-based), kmer_counts, vaf
    """
    # Generate a unique reference sequence
    ref_seq = make_unique_target_sequence(target_length, kmer_length, rng)

    # Choose a variant position away from the edges
    margin = kmer_length + 5
    if margin >= len(ref_seq) // 2:
        margin = max(5, len(ref_seq) // 4)
    pos = rng.randint(margin, len(ref_seq) - margin - max(indel_len, dup_len) - 2)

    # Inject the variant
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

    # Simulate k-mer counts
    kmer_counts = simulate_kmer_counts(
        ref_seq, alt_seq, kmer_length, coverage, vaf, rng
    )

    return {
        "target_id":    target_id,
        "variant_type": variant_type,
        "ref_seq":      ref_seq,
        "alt_seq":      alt_seq,
        "ref_allele":   ref_allele,
        "alt_allele":   alt_allele,
        "pos":          pos + 1,       # 1-based
        "chrom":        f"sim{target_id}",
        "kmer_counts":  kmer_counts,
        "vaf":          vaf,
        "coverage":     coverage,
        "kmer_length":  kmer_length,
    }


# ---------------------------------------------------------------------------
# Main generation logic
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
    indel_min_len: int,
    indel_max_len: int,
    dup_len: int,
    seed: int,
) -> None:
    """Generate a complete simulated benchmark dataset."""
    output_dir.mkdir(parents=True, exist_ok=True)
    targets_dir = output_dir / "targets"
    targets_dir.mkdir(exist_ok=True)

    rng = random.Random(seed)

    # Build variant list
    variant_specs: list[tuple[str, float]] = []

    total_real = n_snvs + n_insertions + n_deletions + n_itds

    if vaf_values:
        # Cycle through explicit VAF values
        vaf_cycle = vaf_values * (total_real // len(vaf_values) + 1)
    else:
        vaf_cycle = None

    def next_vaf(i: int) -> float:
        if vaf_cycle:
            return vaf_cycle[i % len(vaf_cycle)]
        # Log-uniform sampling between vaf_min and vaf_max
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
        variant_specs.append(("Substitution", 0.0))  # absent = VAF 0
        idx += 1

    rng.shuffle(variant_specs)

    # Generate each target
    all_kmer_counts: dict[str, int] = {}
    ground_truth_rows: list[dict] = []
    fasta_records: list[tuple[str, str]] = []

    print(f"Generating {len(variant_specs)} targets in {output_dir}...")

    for i, (vtype, vaf) in enumerate(variant_specs):
        target_id = f"{i+1:04d}"

        # Vary INDEL length
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
                indel_len=indel_len,
                dup_len=dup,
            )
        except Exception as e:
            print(f"  Warning: could not generate target {target_id}: {e}", file=sys.stderr)
            continue

        # Write FASTA
        fasta_path = targets_dir / f"target_{target_id}.fa"
        with open(fasta_path, "w") as f:
            f.write(f">{t['target_id']}\n{t['ref_seq']}\n")
        fasta_records.append((t["target_id"], t["ref_seq"]))

        # Merge k-mer counts (max count per k-mer across targets)
        for kmer, count in t["kmer_counts"].items():
            all_kmer_counts[kmer] = max(all_kmer_counts.get(kmer, 0), count)

        # Ground truth entry
        ground_truth_rows.append({
            "chrom":        t["chrom"],
            "pos":          t["pos"],
            "ref_allele":   t["ref_allele"],
            "alt_allele":   t["alt_allele"],
            "variant_type": vtype,
            "true_vaf":     vaf,
            "category":     vtype,
        })

    # Write mock database TSV
    mock_db_path = output_dir / "mock_db.tsv"
    print(f"Writing {len(all_kmer_counts)} k-mers to mock_db.tsv...")
    with open(mock_db_path, "w") as f:
        f.write("kmer\tcount\n")
        for kmer, count in sorted(all_kmer_counts.items()):
            f.write(f"{kmer}\t{count}\n")

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

    # Write manifest
    n_present = sum(1 for r in ground_truth_rows if r["true_vaf"] > 0)
    manifest = {
        "generated_by": "generate_simulated.py",
        "seed": seed,
        "kmer_length": kmer_length,
        "target_length": target_length,
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
        "total_unique_kmers": len(all_kmer_counts),
        "files": {
            "targets_dir": str(targets_dir),
            "all_targets_fasta": str(all_targets_path),
            "mock_db": str(mock_db_path),
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
    print(f"  Unique k-mers in mock DB: {len(all_kmer_counts)}")
    print(f"  Files:")
    for name, path in manifest["files"].items():
        print(f"    {name}: {path}")
    print()


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

PRESETS = {
    "snv_indel": dict(
        description="Standard SNV + INDEL accuracy dataset",
        n_snvs=30,
        n_insertions=15,
        n_deletions=15,
        n_itds=5,
        n_absent=10,
        vaf_min=0.001,
        vaf_max=0.5,
        coverage=3000,
        kmer_length=31,
        target_length=150,
        indel_min_len=1,
        indel_max_len=10,
    ),
    "large_indels": dict(
        description="Large INDEL sensitivity dataset (5–50 bp)",
        n_snvs=0,
        n_insertions=15,
        n_deletions=15,
        n_itds=5,
        n_absent=5,
        vaf_min=0.01,
        vaf_max=0.2,
        coverage=3000,
        kmer_length=21,
        target_length=200,
        indel_min_len=5,
        indel_max_len=50,
    ),
    "ultra_low_vaf": dict(
        description="Ultra-low VAF dataset (0.01%–0.1%)",
        n_snvs=30,
        n_insertions=10,
        n_deletions=10,
        n_itds=0,
        n_absent=10,
        vaf_min=0.0001,
        vaf_max=0.001,
        coverage=5000,
        kmer_length=31,
        target_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
    "dilution_series": dict(
        description="Fixed VAF levels simulating dilution series",
        n_snvs=10,
        n_insertions=5,
        n_deletions=5,
        n_itds=0,
        n_absent=5,
        vaf_min=None,
        vaf_max=None,
        vaf_values=[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10],
        coverage=3000,
        kmer_length=31,
        target_length=150,
        indel_min_len=1,
        indel_max_len=7,
    ),
}


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
        default=pathlib.Path("data/simulated"),
        help="Base directory for --preset all. Each preset is placed in a subdirectory. "
             "[default: data/simulated/]",
    )
    p.add_argument(
        "--preset",
        choices=list(PRESETS.keys()) + ["all"],
        help="Use a predefined dataset configuration. "
             "Use 'all' to generate all presets under --base-dir.",
    )

    # Variant counts
    p.add_argument("--n-snvs",        type=int, default=30, help="Number of SNVs [default: 30]")
    p.add_argument("--n-insertions",  type=int, default=15, help="Number of insertions [default: 15]")
    p.add_argument("--n-deletions",   type=int, default=15, help="Number of deletions [default: 15]")
    p.add_argument("--n-itds",        type=int, default=5,  help="Number of ITDs [default: 5]")
    p.add_argument("--n-absent",      type=int, default=10,
                   help="Number of absent variants (specificity controls) [default: 10]")

    # VAF settings
    p.add_argument("--vaf-min",   type=float, default=0.001,
                   help="Minimum VAF for log-uniform sampling [default: 0.001]")
    p.add_argument("--vaf-max",   type=float, default=0.5,
                   help="Maximum VAF for log-uniform sampling [default: 0.5]")
    p.add_argument("--vaf-values", type=str,
                   help="Explicit comma-separated VAF values (overrides --vaf-min/max). "
                        "Cycles through for all variants. Example: '0.001,0.01,0.1'")

    # Sequencing parameters
    p.add_argument("--coverage",      type=int, default=3000,
                   help="Simulated sequencing depth (reads, before dedup) [default: 3000]")
    p.add_argument("--kmer-length",   type=int, default=31,
                   help="K-mer length [default: 31]")
    p.add_argument("--target-length", type=int, default=150,
                   help="Target sequence length (bp) [default: 150]")

    # INDEL settings
    p.add_argument("--indel-min-len", type=int, default=1,
                   help="Minimum INDEL length (bp) [default: 1]")
    p.add_argument("--indel-max-len", type=int, default=10,
                   help="Maximum INDEL length (bp) [default: 10]")
    p.add_argument("--dup-len",       type=int, default=10,
                   help="ITD duplication length (bp) [default: 10]")

    p.add_argument("--seed", type=int, default=42, help="Random seed [default: 42]")

    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # Parse explicit VAF values if provided
    vaf_values: Optional[list[float]] = None
    if args.vaf_values:
        vaf_values = [float(v.strip()) for v in args.vaf_values.split(",")]

    if args.preset == "all":
        base = args.base_dir
        for name, cfg in PRESETS.items():
            print(f"=== Generating preset: {name} ===")
            print(f"    {cfg['description']}")
            out = base / name
            generate_dataset(
                output_dir=out,
                n_snvs=cfg["n_snvs"],
                n_insertions=cfg["n_insertions"],
                n_deletions=cfg["n_deletions"],
                n_itds=cfg["n_itds"],
                n_absent=cfg["n_absent"],
                vaf_min=cfg.get("vaf_min") or args.vaf_min,
                vaf_max=cfg.get("vaf_max") or args.vaf_max,
                vaf_values=cfg.get("vaf_values") or vaf_values,
                coverage=cfg["coverage"],
                kmer_length=cfg["kmer_length"],
                target_length=cfg["target_length"],
                indel_min_len=cfg["indel_min_len"],
                indel_max_len=cfg["indel_max_len"],
                dup_len=cfg.get("dup_len", args.dup_len),
                seed=args.seed,
            )
        return 0

    elif args.preset and args.preset in PRESETS:
        cfg = PRESETS[args.preset]
        out = args.output_dir or (args.base_dir / args.preset)
        generate_dataset(
            output_dir=out,
            n_snvs=cfg["n_snvs"],
            n_insertions=cfg["n_insertions"],
            n_deletions=cfg["n_deletions"],
            n_itds=cfg["n_itds"],
            n_absent=cfg["n_absent"],
            vaf_min=cfg.get("vaf_min") or args.vaf_min,
            vaf_max=cfg.get("vaf_max") or args.vaf_max,
            vaf_values=cfg.get("vaf_values") or vaf_values,
            coverage=cfg["coverage"],
            kmer_length=cfg["kmer_length"],
            target_length=cfg["target_length"],
            indel_min_len=cfg["indel_min_len"],
            indel_max_len=cfg["indel_max_len"],
            dup_len=cfg.get("dup_len", args.dup_len),
            seed=args.seed,
        )
        return 0

    else:
        # Custom configuration
        if args.output_dir is None:
            print("ERROR: --output-dir is required (or use --preset)", file=sys.stderr)
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
            indel_min_len=args.indel_min_len,
            indel_max_len=args.indel_max_len,
            dup_len=args.dup_len,
            seed=args.seed,
        )
        return 0


if __name__ == "__main__":
    sys.exit(main())
