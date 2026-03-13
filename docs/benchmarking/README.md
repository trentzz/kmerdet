# kmerdet Benchmarking Framework

This directory contains everything needed to evaluate kmerdet's accuracy,
performance, and scalability: configuration files, runnable scripts, dataset
preparation tools, and templates for writing up comparisons.

---

## What We Benchmark

### 1. Accuracy (Detection Performance)

How well does kmerdet find real variants and avoid false positives?

| Metric | Definition | Baseline (km/kam thesis) |
|--------|-----------|--------------------------|
| SNV sensitivity | TP / (TP + FN) for SNVs | 77% |
| INDEL sensitivity (overall) | TP / (TP + FN) for INDELs | 38% |
| INDEL sensitivity (>7 bp) | TP / (TP + FN) for large INDELs | ~3% |
| Precision (PPV) | TP / (TP + FP) | Not reported in thesis |
| F1-score | Harmonic mean of sensitivity and precision | — |
| VAF detection threshold | Lowest reliably detectable VAF | ~0.1% |

Accuracy is broken down by:
- **Variant type**: Substitution, Insertion, Deletion, ITD, Complex
- **VAF bin**: `[0, 0.1%)`, `[0.1%, 1%)`, `[1%, 5%)`, `[5%, 100%]`
- **kmerdet vs km**: McNemar's test for head-to-head comparison on matched samples

### 2. Performance (Speed and Memory)

How fast is kmerdet compared to km/kmtools, and does it scale well?

| Metric | Baseline (thesis) | Goal |
|--------|------------------|------|
| End-to-end runtime | ~6 min (50 targets, 4 threads) | ≤ 4 min |
| km detect step | ~2 min | < 1 min |
| Peak memory | Not measured | < 4 GB |

Performance is characterized by:
- Runtime vs number of targets (linear scaling expected)
- Runtime vs number of threads (parallelism efficiency)
- Peak RSS memory vs database size
- Throughput: targets per second

### 3. Scalability

Does kmerdet remain fast as the panel grows?

| Panel size | Expected behavior |
|-----------|-------------------|
| 1 target | Sub-second |
| 10 targets | ~10 seconds |
| 50 targets | ≤ 4 minutes |
| 200 targets | ≤ 15 minutes |
| 1000 targets | ≤ 1 hour |

---

## Directory Structure

```
docs/benchmarking/
├── README.md                        ← You are here
│
├── framework/                       ← Runnable benchmark scripts
│   ├── benchmark_config.toml        ← Dataset/parameter configuration
│   ├── run_benchmarks.sh            ← Main entry point for all benchmarks
│   ├── analyze_results.py           ← Compute metrics from TSV output
│   └── plot_results.py              ← Generate figures from metrics
│
├── datasets/                        ← Dataset documentation and generation
│   ├── README.md                    ← How to prepare benchmark datasets
│   └── generate_simulated.py        ← Generate synthetic test data
│
├── results/                         ← Benchmark output (git-ignored large files)
│   └── README.md                    ← Naming conventions and structure
│
├── analysis/                        ← Writeup templates
│   └── comparison_template.md       ← Template for kmerdet vs km analysis
│
└── ci/                              ← Continuous integration
    └── benchmark_ci.yml             ← GitHub Actions workflow
```

---

## How to Run Benchmarks

### Quick Start (5 minutes, simulated data)

```bash
# 1. Generate simulated test data
python3 docs/benchmarking/datasets/generate_simulated.py \
    --output-dir /tmp/kmerdet-bench \
    --n-snvs 20 --n-indels 10 --coverage 3000

# 2. Build kmerdet
cargo build --release

# 3. Run accuracy benchmarks on simulated data
bash docs/benchmarking/framework/run_benchmarks.sh \
    --dataset simulated \
    --data-dir /tmp/kmerdet-bench \
    --output-dir /tmp/kmerdet-results \
    --quick

# 4. Analyze results
python3 docs/benchmarking/framework/analyze_results.py \
    --results-dir /tmp/kmerdet-results \
    --output-dir /tmp/kmerdet-analysis

# 5. Plot results
python3 docs/benchmarking/framework/plot_results.py \
    --metrics /tmp/kmerdet-analysis/metrics.json \
    --output-dir /tmp/kmerdet-analysis/figures
```

### Full Benchmark Suite (requires real data)

```bash
# Edit benchmark_config.toml to point to your datasets
vim docs/benchmarking/framework/benchmark_config.toml

# Run the full suite
bash docs/benchmarking/framework/run_benchmarks.sh \
    --config docs/benchmarking/framework/benchmark_config.toml \
    --output-dir results/$(date +%Y%m%d)
```

### CI Benchmarks (automated, on PR)

Triggered automatically by the GitHub Actions workflow in `ci/benchmark_ci.yml`.
Results are posted as PR comments with a regression table.

---

## How to Interpret Results

### Sensitivity Tables

The `analyze_results.py` script produces tables like:

```
Overall:  sensitivity=0.77  precision=0.92  F1=0.84  (TP=77, FP=7, FN=23)

Per-type breakdown:
  Substitution:  sensitivity=0.77  present=100  TP=77  FN=23
  Insertion:     sensitivity=0.35  present=50   TP=18  FN=32
  Deletion:      sensitivity=0.42  present=50   TP=21  FN=29
  ITD:           sensitivity=0.20  present=10   TP=2   FN=8

VAF bin breakdown:
  [0, 0.001):    sensitivity=0.15  present=30  TP=5   FN=25
  [0.001, 0.01): sensitivity=0.62  present=40  TP=25  FN=15
  [0.01, 0.1):   sensitivity=0.88  present=50  TP=44  FN=6
  [0.1, 1.0]:    sensitivity=0.98  present=80  TP=78  FN=2
```

**What to look for**:
- Sensitivity for Substitution should be >= 77% (thesis baseline).
- INDEL sensitivity should improve over time as algorithmic fixes land.
- Low VAF bins (`[0.001, 0.01)`) are the most important to track for MRD.
- Precision should stay >= 0.90 to avoid clinical false alarms.

### Performance Tables

```
Timing:
  detect (50 targets, 4 threads):  45.2s  ± 1.1s
  filter:                            2.3s  ± 0.1s
  total pipeline:                   47.5s  ± 1.2s

Memory:
  peak_rss_mb:  312
```

A regression is flagged if runtime increases > 20% or memory increases > 50%
compared to the stored baseline.

### ROC-Like Threshold Sweep

The threshold sweep (`--sweep-vaf`) varies the minimum rVAF cutoff and reports
(sensitivity, specificity) at each point. The resulting curve shows the
sensitivity-specificity tradeoff and helps choose an operating point for clinical
use.

A good operating point for MRD monitoring: sensitivity >= 0.80,
specificity >= 0.95.

---

## Glossary

| Term | Meaning |
|------|---------|
| rVAF | Relative Variant Allele Frequency — kmerdet's estimated VAF from NNLS |
| True VAF | Known spike-in or simulation VAF used as ground truth |
| TP | True positive: variant present in truth and detected by kmerdet |
| FP | False positive: detected by kmerdet but not in truth (or truth VAF = 0) |
| FN | False negative: present in truth but not detected by kmerdet |
| TN | True negative: absent in truth and not detected |
| MCC | Matthews Correlation Coefficient: balanced metric for imbalanced classes |
| F1 | Harmonic mean of precision and sensitivity |
| McNemar | Paired statistical test comparing kmerdet vs km on same samples |

---

## Dependencies

### For benchmark execution
- `kmerdet` binary (Rust, built with `cargo build --release`)
- `jellyfish` >= 2.3 (for real data; not needed for simulated mock data)
- `hyperfine` (optional, for accurate timing; falls back to `/usr/bin/time`)

### For analysis and plotting
- Python >= 3.9
- `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`

Install Python dependencies:
```bash
pip install numpy pandas scipy scikit-learn matplotlib
```

---

## Adding a New Dataset

1. Place your jellyfish database (`.jf`) and target FASTA files in a directory.
2. Create a ground truth TSV (format: see `datasets/README.md`).
3. Add a `[[dataset]]` entry to `framework/benchmark_config.toml`.
4. Run `run_benchmarks.sh --dataset <name>`.

---

## Contributing Results

After running benchmarks, save results with a timestamped name:

```bash
cp -r /tmp/kmerdet-results docs/benchmarking/results/$(date +%Y%m%d)-v0.1.0-simulated/
```

Commit the `metrics.json` and summary TSVs (not raw TSV outputs, which can be
large). The `results/README.md` has naming conventions.
