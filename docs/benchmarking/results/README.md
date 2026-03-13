# Benchmark Results

This directory stores benchmark results from kmerdet runs over time.

Large raw output files (TSV detection results, full per-variant tables) are
git-ignored via `.gitignore`. Only the compact summary files should be
committed.

---

## Naming Convention

Each benchmark run is stored in a timestamped subdirectory:

```
results/
├── README.md                        ← This file
├── baseline_km_thesis.json          ← Reference: km/kam thesis performance
├── YYYYMMDD-vX.Y.Z-datasetname/     ← One directory per run
│   ├── run_metadata.json            ← Date, version, git hash, parameters
│   ├── aggregate_accuracy.json      ← All dataset metrics in one file
│   ├── snv_indel_benchmark.json     ← Per-dataset benchmark JSON
│   ├── large_indels_benchmark.json
│   ├── performance.tsv              ← Timing and memory
│   └── figures/                     ← Plots (optional)
│       ├── sensitivity_by_vaf.png
│       ├── variant_type_breakdown.png
│       └── sensitivity_heatmap.png
└── latest -> YYYYMMDD-vX.Y.Z.../   ← Symlink to most recent run
```

### Directory name format

```
YYYYMMDD-vMAJOR.MINOR.PATCH-DATASETNAME
```

Examples:
- `20240115-v0.1.0-simulated_snv_indel`
- `20240302-v0.2.0-patient_cohort_10`
- `20240415-v0.3.0-full`

### When to create a new result directory

- After each numbered kmerdet release (v0.1.0, v0.2.0, etc.)
- After any algorithmic change that could affect sensitivity or precision
- After adding a new dataset to the benchmark suite

---

## Files to Commit

Commit these (small, human-readable):
- `run_metadata.json`
- `aggregate_accuracy.json`
- `*_benchmark.json` (from kmerdet benchmark subcommand)
- `performance.tsv`
- Figure files (PNG/PDF)

Do NOT commit these (large):
- Raw detection TSV files (`.tsv` with full variant calls)
- Intermediate files from run_benchmarks.sh

---

## Baseline Reference: km/kam Thesis Performance

The file `baseline_km_thesis.json` records the expected performance of the
original km/kam pipeline from the thesis validation study (n=10 patients).
All kmerdet results are compared against this baseline.

```json
{
  "source": "Thesis validation study, n=10 patients, UMI duplex sequencing",
  "cohort": "10 patients, serial blood draws, patient-specific targeted panels",
  "pipeline": "HUMID + jellyfish + kmtools chunk/filter",
  "parameters": {
    "k": 31,
    "count": 2,
    "ratio": 0.00001,
    "flank": 35,
    "threads": 4
  },
  "performance": {
    "overall_sensitivity":   0.60,
    "snv_sensitivity":       0.77,
    "indel_sensitivity_overall": 0.38,
    "indel_sensitivity_large":   0.03,
    "vaf_detection_threshold":   0.001,
    "runtime_seconds_50targets": 360,
    "peak_rss_mb":               null
  },
  "notes": [
    "SNV sensitivity 77% vs 90-95% for alignment-based pipelines",
    "INDEL sensitivity 38% overall; near-zero for INDELs >7 bp",
    "Detection threshold ~0.1% VAF (count=2 at 3000x coverage)",
    "6 min total including HUMID dedup and jellyfish count",
    "2 min for km detect step on 50-target panel with 4 threads"
  ]
}
```

---

## Interpreting Results Over Time

### Tracking regressions

Run the comparison after each major commit:

```bash
python3 docs/benchmarking/framework/analyze_results.py \
    --results-dir docs/benchmarking/results/latest/ \
    --output-dir /tmp/analysis

# Compare to previous run
python3 - <<'EOF'
import json, pathlib

current = json.load(open("/tmp/analysis/metrics.json"))
baseline = json.load(open("docs/benchmarking/results/baseline_km_thesis.json"))

curr_snv = current.get("per_type", [{}])[0].get("sensitivity") or \
           current.get("overall", {}).get("sensitivity", 0)
base_snv = baseline["performance"]["snv_sensitivity"]

delta = curr_snv - base_snv
print(f"SNV sensitivity: {curr_snv:.3f} (baseline: {base_snv:.3f}, delta: {delta:+.3f})")
if delta < -0.05:
    print("WARNING: Sensitivity regression > 5 pp vs baseline!")
elif delta > 0.05:
    print("Improvement: > 5 pp gain over baseline.")
else:
    print("Within expected range.")
EOF
```

### Expected progression

As algorithmic improvements from the roadmap are implemented:

| kmerdet version | Expected SNV sensitivity | Expected INDEL sensitivity | Notes |
|----------------|--------------------------|----------------------------|-------|
| v0.1.0 (Phase 1-3) | ≥ 77% (match km) | ≥ 38% (match km) | Correctness baseline |
| v0.2.0 (adaptive thresholds) | ≥ 82% | ≥ 45% | Adaptive count/ratio |
| v0.3.0 (multi-k) | ≥ 82% | ≥ 65% | Per-target k optimization |
| v0.4.0 (UMI-aware counting) | ≥ 85% | ≥ 70% | Molecular-level evidence |
| v1.0.0 (full pipeline) | ≥ 87% | ≥ 75% | All improvements integrated |

These are targets, not guarantees. The benchmarks will show whether each
feature delivers the expected improvement.

---

## Comparing Across Datasets

When comparing results across datasets (e.g., simulated vs real), note that:

1. **Simulated data is optimistic**: no GC bias, no repetitive regions, no
   off-target k-mers. Real sensitivity on clinical data will be lower.

2. **VAF distribution matters**: a dataset with more ultra-low-VAF variants will
   show lower overall sensitivity than one with mostly high-VAF variants, even if
   the algorithm is identical.

3. **Variant type mix matters**: a dataset with many large INDELs will show lower
   overall sensitivity due to the k-mer length constraint.

4. **Sequencing depth matters**: the simulated datasets use 3000x coverage by
   default. Clinical samples may have higher or lower depth.

Always report the dataset composition alongside sensitivity numbers.
